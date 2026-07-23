(ns evophy.qm.matrix
  "Second, independent oracle for the 1D time-independent Schrodinger equation:
  finite-difference Hamiltonian + dense symmetric diagonalization (cyclic
  Jacobi rotations), instead of Numerov + shooting.

  This exists to cross-validate `evophy.qm.schrodinger`, not to replace it: a
  bug in the only oracle would silently poison the fitness layer and the
  search built on top of it, so a second method built on entirely different
  numerics (no integration, no bisection, no shooting) is a check the first
  method can't check itself.

  Method: discretize -(hbar^2/2m) psi''(x) + V(x) psi(x) = E psi(x) on the
  interior grid points (Dirichlet psi=0 at both boundaries, matching the
  shooting solver's convention) with the standard 3-point second-derivative
  stencil. This gives a dense (in practice tridiagonal, but we diagonalize it
  as a general symmetric matrix) eigenvalue problem H psi = E psi whose
  eigenpairs are exactly the bound-state energies/wavefunctions on that grid.

  Pure Clojure, no linear-algebra dependency (same constraint as
  `schrodinger.clj`): the eigensolver is a hand-rolled classical Jacobi
  rotation sweep. That's O(n^3) per sweep, so this module uses much smaller
  grids than the shooting solver's default (a few hundred points, not
  thousands) -- plenty for cross-validation, not for production-scale search."
  (:require [evophy.qm.schrodinger :as s]))

;; ---------------------------------------------------------------------------
;; Hamiltonian matrix (finite-difference, interior points only)
;; ---------------------------------------------------------------------------

(defn hamiltonian
  "Dense symmetric Hamiltonian for the interior grid points of `params`
  (indices 1..n-2 of `(s/grid params)`; the boundary points are pinned to
  psi=0 and excluded). Row/column i of the returned matrix corresponds to
  grid point i+1.

  Diagonal: 2k + V(x_i); off-diagonal: -k, where k = hbar^2/(2 m h^2) is the
  finite-difference kinetic coefficient. Returns a vector of double-arrays."
  [{:keys [hbar m V] :as params}]
  (let [xs (s/grid params)
        h  (s/step-size params)
        k  (/ (* hbar hbar) (* 2.0 m h h))
        n  (count xs)
        ni (- n 2)]
    (vec (for [i (range ni)]
           (let [row (double-array ni)
                 xi  (nth xs (inc i))]
             (aset row i (+ (* 2.0 k) (double (V xi))))
             (when (pos? i) (aset row (dec i) (- k)))
             (when (< i (dec ni)) (aset row (inc i) (- k)))
             row)))))

;; ---------------------------------------------------------------------------
;; Cyclic Jacobi eigenvalue algorithm (dense symmetric matrices)
;; ---------------------------------------------------------------------------

(defn- copy-matrix
  ^objects [rows]
  (into-array (map #(aclone ^doubles %) rows)))

(defn- identity-matrix
  ^objects [^long n]
  (into-array
    (for [i (range n)]
      (let [row (double-array n)]
        (aset row i 1.0)
        row))))

(defn- mget ^double [^objects m ^long i ^long j]
  (aget ^doubles (aget m i) j))

(defn- mset! [^objects m ^long i ^long j ^double v]
  (aset ^doubles (aget m i) j v))

(defn- frobenius-norm
  "sqrt(sum of squares of all entries), used to scale the convergence check."
  ^double [^objects a ^long n]
  (loop [i 0 s 0.0]
    (if (< i n)
      (recur (inc i)
             (+ s (loop [j 0 s2 0.0]
                    (if (< j n)
                      (recur (inc j) (+ s2 (let [v (mget a i j)] (* v v))))
                      s2))))
      (Math/sqrt s))))

(defn- off-diagonal-norm
  ^double [^objects a ^long n]
  (loop [i 0 s 0.0]
    (if (< i n)
      (recur (inc i)
             (+ s (loop [j (inc i) s2 0.0]
                    (if (< j n)
                      (recur (inc j) (+ s2 (let [v (mget a i j)] (* v v))))
                      s2))))
      (Math/sqrt (* 2.0 s)))))

(defn- jacobi-rotate!
  "Zero A[p][q] (and A[q][p]) in-place via a single Jacobi rotation, updating
  the rest of A and accumulating the rotation into V (whose columns converge
  to the eigenvectors)."
  [^objects a ^objects v n p q]
  (let [n   (long n)
        p   (long p)
        q   (long q)
        apq (mget a p q)]
    (when-not (zero? apq)
      (let [app   (mget a p p)
            aqq   (mget a q q)
            theta (/ (- aqq app) (* 2.0 apq))
            sign  (if (neg? theta) -1.0 1.0)
            t     (/ sign (+ (Math/abs theta) (Math/sqrt (+ 1.0 (* theta theta)))))
            c     (/ 1.0 (Math/sqrt (+ 1.0 (* t t))))
            s     (* t c)
            tau   (/ s (+ 1.0 c))]
        (mset! a p p (- app (* t apq)))
        (mset! a q q (+ aqq (* t apq)))
        (mset! a p q 0.0)
        (mset! a q p 0.0)
        (dotimes [i n]
          (when (and (not= i p) (not= i q))
            (let [aip (mget a i p)
                  aiq (mget a i q)
                  nip (- aip (* s (+ aiq (* tau aip))))
                  niq (+ aiq (* s (- aip (* tau aiq))))]
              (mset! a i p nip) (mset! a p i nip)
              (mset! a i q niq) (mset! a q i niq))))
        (dotimes [i n]
          (let [vip (mget v i p)
                viq (mget v i q)]
            (mset! v i p (- vip (* s (+ viq (* tau vip)))))
            (mset! v i q (+ viq (* s (- vip (* tau viq)))))))))))

(defn jacobi-eigen
  "Diagonalize dense symmetric matrix `rows` (n x n, vector of double-arrays)
  via cyclic Jacobi sweeps. Returns {:eigenvalues (double-array n, unsorted)
  :eigenvectors (n x n objects, column j is the unit eigenvector for
  eigenvalues[j]) :sweeps (int)}.

  Converges quadratically once off-diagonal mass is small; `tol` is relative
  to the matrix's Frobenius norm so it's scale-invariant."
  [rows n & {:keys [max-sweeps tol] :or {max-sweeps 100 tol 1e-13}}]
  (let [n     (long n)
        a     (copy-matrix rows)
        v     (identity-matrix n)
        scale (max (frobenius-norm a n) 1e-300)]
    (loop [sweep 0]
      (if (or (>= sweep max-sweeps) (< (/ (off-diagonal-norm a n) scale) tol))
        (let [d (double-array n)]
          (dotimes [i n] (aset d i (mget a i i)))
          {:eigenvalues d :eigenvectors v :sweeps sweep})
        (do
          (dotimes [p n]
            (dotimes [dq (- n p 1)]
              (jacobi-rotate! a v n p (+ p dq 1))))
          (recur (inc sweep)))))))

;; ---------------------------------------------------------------------------
;; Bound-state solve (finite-difference + diagonalize)
;; ---------------------------------------------------------------------------

(defn- count-nodes-physical
  "Sign changes restricted to the classically-allowed region V(x) <= E, same
  convention as `schrodinger.clj`'s node counter. Needed because a Jacobi
  eigenvector's classically-forbidden tail decays toward true zero, so
  roundoff noise (not physics) dominates its sign there -- counting flips in
  that noise floor would wildly overcount nodes on anything but a hard wall."
  ^long [{:keys [V]} xs psi ^double E]
  (let [n (count psi)]
    (loop [i 2 prev (double (nth psi 1)) nodes 0]
      (if (< i n)
        (let [allowed? (<= (double (V (nth xs i))) E)
              v (double (nth psi i))]
          (recur (inc i)
                 (if (zero? v) prev v)
                 (if (and allowed? (not (zero? v)) (neg? (* prev v)))
                   (inc nodes) nodes)))
        nodes))))

(defn solve
  "Find the first `n-levels` bound states of `params` by finite-difference +
  dense diagonalization. Returns a vector of {:level :energy :psi :nodes}
  maps ordered by energy, same shape as `evophy.qm.schrodinger/solve`.

  `psi` is padded with the Dirichlet zeros at both ends and unit-normalized
  under the trapezoidal L2 inner product (matching the oracle's convention):
  a Jacobi eigenvector is Euclidean-unit (sum v_i^2 = 1), and since the
  padded boundary points contribute 0, L2-norm^2 = h * (Euclidean norm)^2 = h,
  so dividing by sqrt(h) makes it L2-unit."
  [params & {:keys [n-levels] :or {n-levels 4}}]
  (let [h    (s/step-size params)
        xs   (s/grid params)
        n    (count xs)
        ni   (- n 2)
        H    (hamiltonian params)
        {:keys [eigenvalues eigenvectors]} (jacobi-eigen H ni)
        order (sort-by #(aget ^doubles eigenvalues %) (range ni))
        scale (/ 1.0 (Math/sqrt h))]
    (vec (for [level (range (min n-levels ni))
               :let [col (nth order level)
                     energy (aget ^doubles eigenvalues col)
                     psi (vec (concat [0.0]
                                      (for [i (range ni)]
                                        (* scale (mget ^objects eigenvectors i col)))
                                      [0.0]))]]
           {:level level :energy energy :psi psi
            :nodes (count-nodes-physical params xs psi energy)}))))

(defn benchmark-report
  "Solve a `schrodinger.clj` benchmark map with the matrix method and compare
  to its closed form. Mirrors `schrodinger/benchmark-report`'s shape so the
  two oracles can be diffed directly."
  [{:keys [params exact search]} & {:keys [n-levels]}]
  (let [n-levels (or n-levels (:n-levels search) 4)
        levels (solve params :n-levels n-levels)]
    (mapv (fn [{:keys [level energy]}]
            (let [e (exact level)]
              {:level level
               :numeric energy
               :exact e
               :rel-error (/ (Math/abs (- energy e)) (Math/abs e))}))
          levels)))

(defn spectrum
  "Solve a benchmark map for its bound-state energies (no closed form needed)."
  [{:keys [params search]} & {:keys [n-levels]}]
  (solve params :n-levels (or n-levels (:n-levels search) 4)))

;; ---------------------------------------------------------------------------
;; Demo: cross-check against the closed forms and against the shooting solver
;; ---------------------------------------------------------------------------

(def ^:private matrix-n
  "Grid size for the demo/benchmarks here. Dense diagonalization is O(n^3),
  so this is far smaller than the shooting solver's default (2001) -- plenty
  of accuracy for cross-validation, not meant for production-scale search."
  301)

(defn -main
  [& _]
  (println "Infinite square well (L=1, hbar=m=1):  E_n = n^2 pi^2 / 2")
  (doseq [row (benchmark-report (s/infinite-well {:L 1.0 :n matrix-n}))]
    (println (format "  level %d: numeric=%.6f exact=%.6f rel-err=%.2e"
                     (:level row) (:numeric row) (:exact row) (:rel-error row))))
  (println)
  (println "Harmonic oscillator (omega=1, hbar=m=1):  E_n = (n + 1/2)")
  (doseq [row (benchmark-report (s/harmonic {:omega 1.0 :n matrix-n}))]
    (println (format "  level %d: numeric=%.6f exact=%.6f rel-err=%.2e"
                     (:level row) (:numeric row) (:exact row) (:rel-error row))))
  (println)
  (println "Finite square well (V0=20, a=1):  transcendental closed form")
  (doseq [row (benchmark-report (s/finite-well {:V0 20.0 :a 1.0 :n matrix-n}))]
    (println (format "  level %d: numeric=%.6f exact=%.6f rel-err=%.2e"
                     (:level row) (:numeric row) (:exact row) (:rel-error row))))
  (println)
  (println "Woods-Saxon (V0=50, R=3, a=0.6): matrix vs shooting (no closed form)")
  (let [bench (s/woods-saxon {:n matrix-n})
        mat   (spectrum bench :n-levels 6)
        shoot (s/spectrum bench :n-levels 6)]
    (doseq [[m sh] (map vector mat shoot)]
      (println (format "  level %d: matrix=%.6f shooting=%.6f diff=%.2e"
                       (:level m) (:energy m) (:energy sh)
                       (Math/abs (- (:energy m) (:energy sh))))))))
