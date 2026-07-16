(ns evophy.qm.schrodinger
  "Ground truth for the 1D time-independent Schrodinger equation

       -(hbar^2 / 2m) psi''(x) + V(x) psi(x) = E psi(x).

  This is the QM analogue of evophy's symplectic integrator for Kepler:
  a numerical oracle that produces checkable answers (bound-state energies
  and wavefunctions) for a given potential V(x). The GP/analytical search
  layered on top scores trial energies / wavefunctions / strategies against
  what this solver returns.

  Method: Numerov integration + shooting. Pure Clojure, no linear-algebra
  dependency. We rewrite the ODE as

       psi''(x) + f(x) psi(x) = 0,   f(x) = (2m/hbar^2)(E - V(x))

  and use the Numerov three-term recurrence, which is O(h^4) accurate for
  exactly this f*psi form. For a trial energy E we integrate across the grid
  and read off the boundary mismatch; bound-state energies are the E values
  where the wavefunction satisfies the far boundary condition (psi -> 0).
  We locate them by counting interior nodes (which increases by one at each
  eigenvalue) and bisecting.

  Units: natural units with hbar = 1, m = 1 by default. Every entry point
  takes an explicit params map so nuclear-scale constants can be supplied
  later without touching the solver."
  (:require [clojure.string :as str]))

;; ---------------------------------------------------------------------------
;; Grid + parameters
;; ---------------------------------------------------------------------------

(def default-params
  "Default physical + numerical parameters. `V` is a fn of x returning the
  potential energy. Override any key per problem."
  {:hbar 1.0
   :m    1.0
   :x-min -10.0
   :x-max  10.0
   :n      2001            ;; grid points (odd => symmetric about origin)
   :V     (fn [_] 0.0)})

(defn grid
  "Vector of x sample points for the given params."
  [{:keys [x-min x-max n]}]
  (let [h (/ (- x-max x-min) (double (dec n)))]
    (mapv (fn [i] (+ x-min (* i h))) (range n))))

(defn step-size
  [{:keys [x-min x-max n]}]
  (/ (- x-max x-min) (double (dec n))))

;; ---------------------------------------------------------------------------
;; Numerov integration for a trial energy
;; ---------------------------------------------------------------------------

(defn- f-array
  "f(x) = (2m/hbar^2)(E - V(x)) sampled on the grid."
  ^doubles [{:keys [hbar m V]} xs ^double E]
  (let [c (/ (* 2.0 m) (* hbar hbar))
        n (count xs)
        out (double-array n)]
    (dotimes [i n]
      (aset out i (* c (- E (double (V (nth xs i)))))))
    out))

(defn integrate
  "Numerov-integrate psi across the grid for trial energy E.

  Left boundary condition: psi[0] = 0, psi[1] = seed (a small positive
  value; overall scale is arbitrary for an eigenvalue search). Returns a
  map with the full :psi double-array and :end (psi at the last grid point),
  the quantity whose sign/zero-crossings drive the shooting search."
  [params xs ^double E]
  (let [h  (step-size params)
        h2 (* h h)
        f  (f-array params xs E)
        n  (count xs)
        psi (double-array n)
        seed 1.0e-4
        c12 (/ h2 12.0)]
    (aset psi 0 0.0)
    (aset psi 1 seed)
    (loop [i 1]
      (if (< i (dec n))
        (let [fm (aget f (dec i))
              fi (aget f i)
              fp (aget f (inc i))
              w-m (* (aget psi (dec i)) (+ 1.0 (* c12 fm)))
              w-i (* (aget psi i)       (- 1.0 (* 5.0 c12 fi)))
              next (/ (- (* 2.0 w-i) w-m)
                      (+ 1.0 (* c12 fp)))]
          (aset psi (inc i) next)
          (recur (inc i)))
        {:psi psi :end (aget psi (dec n))}))))

(defn- count-nodes-full
  "Total sign changes of a left-shot psi across the whole grid, including the
  final index. This drives the eigenvalue bisection: as trial E rises past a
  bound-state energy, one extra sign change appears — at the wall for a
  hard-wall potential, or in the runaway forbidden tail for a decaying state.
  In both cases the transition coincides with the eigenvalue, so bisecting on
  this count converges to it. (The tail crossing is spurious as *physics* but
  is exactly the shooting signal, hence used for bracketing only.)"
  ^long [^doubles psi]
  (let [n (alength psi)]
    (loop [i 2 prev (aget psi 1) nodes 0]
      (if (< i n)
        (let [v (aget psi i)]
          (recur (inc i)
                 (if (zero? v) prev v)
                 (if (and (not (zero? v)) (neg? (* prev v)))
                   (inc nodes) nodes)))
        nodes))))

(defn- count-nodes-physical
  "Physical node count: sign changes only in the classically-allowed region
  V(x) <= E, ignoring the runaway forbidden tail. Equals the quantum number
  of the bound state; used for reporting/validation, not for bracketing."
  ^long [{:keys [V]} xs ^doubles psi ^double E]
  (let [n (alength psi)]
    (loop [i 2 prev (aget psi 1) nodes 0]
      (if (< i n)
        (let [allowed? (<= (double (V (nth xs i))) E)
              v (aget psi i)]
          (recur (inc i)
                 (if (zero? v) prev v)
                 (if (and allowed? (not (zero? v)) (neg? (* prev v)))
                   (inc nodes) nodes)))
        nodes))))

;; ---------------------------------------------------------------------------
;; Eigenvalue search (shooting + bisection)
;; ---------------------------------------------------------------------------

(defn- normalize
  "Trapezoidal-normalize psi to unit L2 norm on the grid, returned as a vector."
  [params ^doubles psi]
  (let [h (step-size params)
        n (alength psi)
        norm2 (loop [i 0 s 0.0]
                (if (< i n)
                  (let [v (aget psi i)
                        w (if (or (zero? i) (= i (dec n))) 0.5 1.0)]
                    (recur (inc i) (+ s (* w v v h))))
                  s))
        scale (/ 1.0 (Math/sqrt (max norm2 1e-300)))]
    (mapv (fn [i] (* scale (aget psi i))) (range n))))

(defn- turning-index
  "Largest grid index still classically allowed (V(x) <= E): the outer turning
  point. Clamped to a usable interior range."
  ^long [{:keys [V]} xs ^double E]
  (let [n (count xs)]
    (loop [i (- n 2) found 1]
      (if (< i 1)
        found
        (if (<= (double (V (nth xs i))) E)
          (max 1 (min i (- n 2)))
          (recur (dec i) found))))))

(defn- integrate-inward
  "Numerov-integrate from the right boundary (psi[n-1]=0, psi[n-2]=seed,
  decaying) leftward down to index `m`. Returns the double-array with indices
  [m .. n-1] filled. Used to give the classically-forbidden tail a stable,
  decaying shape instead of the runaway a left-only shot produces."
  ^doubles [params xs ^double E ^long m]
  (let [h2 (let [h (step-size params)] (* h h))
        f  (f-array params xs E)
        n  (count xs)
        psi (double-array n)
        c12 (/ h2 12.0)]
    (aset psi (dec n) 0.0)
    (aset psi (- n 2) 1.0e-4)
    (loop [i (- n 2)]
      (if (> i m)
        (let [fi  (aget f i)
              fp  (aget f (inc i))
              fm  (aget f (dec i))
              w-i (* (aget psi i)       (- 1.0 (* 5.0 c12 fi)))
              w-p (* (aget psi (inc i)) (+ 1.0 (* c12 fp)))
              prev (/ (- (* 2.0 w-i) w-p)
                      (+ 1.0 (* c12 fm)))]
          (aset psi (dec i) prev)
          (recur (dec i)))
        psi))))

(defn- build-wavefunction
  "Construct a clean, unit-normalized wavefunction at (near-)eigenvalue E.

  Shooting only from the left gives a wavefunction that runs away in the
  classically-forbidden tail. Instead we glue the left-shot solution (valid up
  to the outer turning point m) to a right-shot decaying solution (valid in the
  tail), amplitude-matched at m. Since E is already resolved to ~1e-10, the
  slope mismatch at the join is negligible, and the result is normalizable."
  [params xs ^double E]
  (let [n (count xs)
        m (turning-index params xs E)
        psi-l (:psi (integrate params xs E))]
    (if (>= m (- n 2))
      ;; Allowed all the way to the wall (e.g. infinite well): no runaway tail.
      (normalize params psi-l)
      (let [psi-r (integrate-inward params xs E m)
            lm (aget ^doubles psi-l m)
            rm (aget ^doubles psi-r m)
            scale (if (zero? rm) 0.0 (/ lm rm))
            glued (double-array n)]
        (dotimes [i n]
          (aset glued i (if (<= i m)
                          (aget ^doubles psi-l i)
                          (* scale (aget ^doubles psi-r i)))))
        (normalize params glued)))))

(defn solve-level
  "Find the bound state with exactly `level` nodes (level 0 = ground state)
  within [E-lo, E-hi] by bisecting on the classically-allowed node count,
  which increases by one at each eigenvalue. Robust for hard walls and
  decaying tails alike. The wavefunction is rebuilt two-sided (see
  `build-wavefunction`) so its tail decays instead of running away.

  Returns {:energy :psi :nodes} or nil if the level is not bracketed."
  [params xs level E-lo E-hi
   & {:keys [tol max-iter] :or {tol 1e-12 max-iter 200}}]
  (let [nodes-of (fn [E] (count-nodes-full (:psi (integrate params xs E))))]
    (when (and (<= (nodes-of E-lo) level) (> (nodes-of E-hi) level))
      (loop [lo E-lo hi E-hi i 0]
        (let [mid (* 0.5 (+ lo hi))
              nd (nodes-of mid)]
          (if (or (> i max-iter) (< (- hi lo) tol))
            {:energy mid
             :psi (build-wavefunction params xs mid)
             :nodes (count-nodes-physical params xs (:psi (integrate params xs mid)) mid)}
            ;; Levels <= target sit below `mid`; raise lo when mid still has
            ;; too few nodes, else lower hi.
            (if (<= nd level)
              (recur mid hi (inc i))
              (recur lo mid (inc i)))))))))

(defn solve
  "Find the first `n-levels` bound states of the potential in `params`.

  Scans trial energies from `E-lo` to `E-hi` in `scan-steps` increments,
  brackets each eigenvalue by watching the interior node count increase,
  and refines each with `solve-level`. Returns a vector of
  {:level :energy :psi :nodes} maps ordered by energy."
  [params & {:keys [n-levels E-lo E-hi scan-steps]
             :or {n-levels 4 E-lo 0.0 E-hi 50.0 scan-steps 4000}}]
  (let [xs (grid params)
        dE (/ (- E-hi E-lo) (double scan-steps))
        node-at (fn [E] (count-nodes-full (:psi (integrate params xs E))))]
    (loop [i 0
           prev-nodes (node-at E-lo)
           found []]
      (if (or (>= (count found) n-levels) (> i scan-steps))
        found
        (let [E (+ E-lo (* (inc i) dE))
              nd (node-at E)]
          (if (> nd prev-nodes)
            ;; A new node appeared between (i-1) and i: eigenvalue for
            ;; `prev-nodes` is bracketed in [E-dE, E].
            (let [lvl prev-nodes
                  hit (solve-level params xs lvl (- E dE) E)]
              (recur (inc i) nd (if hit (conj found (assoc hit :level lvl)) found)))
            (recur (inc i) nd found)))))))

;; ---------------------------------------------------------------------------
;; Benchmark potentials with closed-form spectra (validation targets)
;; ---------------------------------------------------------------------------

(defn infinite-well
  "Infinite square well of width L on [0, L]. Represent the hard walls with a
  grid spanning exactly [0, L] and Dirichlet psi(0)=psi(L)=0 (automatic here).

  Closed form: E_n = n^2 * pi^2 * hbar^2 / (2 m L^2), n = 1,2,3,..."
  [{:keys [hbar m L n] :or {hbar 1.0 m 1.0 L 1.0 n 2001}}]
  {:params (merge default-params
                  {:hbar hbar :m m :x-min 0.0 :x-max L :n n
                   :V (fn [_] 0.0)})
   :exact  (fn [^long lvl]            ;; lvl 0 => n=1
             (let [nn (inc lvl)]
               (/ (* nn nn Math/PI Math/PI hbar hbar)
                  (* 2.0 m L L))))})

(defn harmonic
  "Harmonic oscillator V(x) = 1/2 m omega^2 x^2.

  Closed form: E_n = hbar * omega * (n + 1/2), n = 0,1,2,..."
  [{:keys [hbar m omega x-max n] :or {hbar 1.0 m 1.0 omega 1.0 x-max 10.0 n 2001}}]
  {:params (merge default-params
                  {:hbar hbar :m m :x-min (- x-max) :x-max x-max :n n
                   :V (fn [x] (* 0.5 m omega omega x x))})
   :exact  (fn [^long lvl]
             (* hbar omega (+ lvl 0.5)))
   :search {:E-lo 0.0 :E-hi (* 12.0 hbar omega) :n-levels 6}})

;; -- transcendental root helpers (for the finite-well closed form) ----------

(defn- bisect-root
  "Bisect for a root of f on [a,b] where f(a) and f(b) straddle zero."
  [f a b tol max-iter]
  (loop [lo a hi b flo (f a) i 0]
    (let [mid (* 0.5 (+ lo hi)) fm (f mid)]
      (if (or (> i max-iter) (< (- hi lo) tol))
        mid
        (if (neg? (* flo fm))
          (recur lo mid flo (inc i))
          (recur mid hi fm (inc i)))))))

(defn- scan-roots
  "All roots of f on [a,b], located by sampling `steps` subintervals and
  bisecting each sign change."
  [f a b steps tol]
  (let [dx (/ (- b a) (double steps))]
    (loop [i 0 prev (f a) roots []]
      (if (>= i steps)
        roots
        (let [x0 (+ a (* i dx))
              x1 (+ a (* (inc i) dx))
              fx1 (f x1)]
          (recur (inc i) fx1
                 (if (and (not (zero? prev)) (neg? (* prev fx1)))
                   (conj roots (bisect-root f x0 x1 tol 200))
                   roots)))))))

(defn finite-well-levels
  "Exact bound-state energies (E<0, ascending) of a 1D finite square well of
  depth V0 and half-width a: V(x) = -V0 for |x|<a, 0 outside; hbar=m=1.

  Inside wavenumber k=sqrt(2(E+V0)), outside decay kappa=sqrt(-2E). Even parity
  states satisfy k*tan(ka)=kappa, odd states k*cot(ka)=-kappa. We use the
  singularity-free forms k sin(ka) - kappa cos(ka) = 0 (even) and
  k cos(ka) + kappa sin(ka) = 0 (odd) and collect all roots in (-V0, 0)."
  [V0 a]
  (let [k   (fn [E] (Math/sqrt (* 2.0 (+ (double E) V0))))
        kap (fn [E] (Math/sqrt (* -2.0 (double E))))
        f-even (fn [E] (- (* (k E) (Math/sin (* (k E) a)))
                          (* (kap E) (Math/cos (* (k E) a)))))
        f-odd  (fn [E] (+ (* (k E) (Math/cos (* (k E) a)))
                          (* (kap E) (Math/sin (* (k E) a)))))
        lo (+ (- V0) 1e-9)
        hi (- 1e-9)]
    (vec (sort (concat (scan-roots f-even lo hi 4000 1e-13)
                       (scan-roots f-odd  lo hi 4000 1e-13))))))

(defn finite-well
  "Finite square well of depth V0 and half-width a on [-x-max, x-max].
  Closed form via `finite-well-levels` (transcendental)."
  [{:keys [V0 a x-max n] :or {V0 20.0 a 1.0 x-max 6.0 n 2001}}]
  (let [levels (finite-well-levels V0 a)]
    {:params (merge default-params
                    {:x-min (- x-max) :x-max x-max :n n
                     :V (fn [x] (if (< (Math/abs (double x)) a) (- V0) 0.0))})
     :exact  (fn [^long lvl] (nth levels lvl))
     :search {:E-lo (- V0) :E-hi -1e-3 :n-levels (count levels)}}))

(defn woods-saxon
  "1D Woods-Saxon well V(x) = -V0 / (1 + exp((|x|-R)/a)) — the standard smooth
  nuclear mean-field shape (symmetric here in |x|). NO closed-form spectrum:
  this is a genuine target, validated only by grid convergence / internal
  consistency, not against an analytic answer."
  [{:keys [V0 R a x-max n] :or {V0 50.0 R 3.0 a 0.6 x-max 12.0 n 2001}}]
  {:params (merge default-params
                  {:x-min (- x-max) :x-max x-max :n n
                   :V (fn [x] (/ (- V0)
                                 (+ 1.0 (Math/exp (/ (- (Math/abs (double x)) R) a)))))})
   :search {:E-lo (- V0) :E-hi -1e-2 :n-levels 8}})

;; ---------------------------------------------------------------------------
;; Convenience: solve a benchmark and report against its closed form
;; ---------------------------------------------------------------------------

(defn benchmark-report
  "Solve a benchmark map and compare to its closed form. Search range/level
  count come from the map's :search hints (overridable via kwargs). Requires
  an :exact fn; use `spectrum` for closed-form-less potentials."
  [{:keys [params exact search]} & {:keys [n-levels E-lo E-hi]}]
  (let [E-lo (or E-lo (:E-lo search) 0.0)
        E-hi (or E-hi (:E-hi search) 200.0)
        n-levels (or n-levels (:n-levels search) 4)
        levels (solve params :n-levels n-levels :E-lo E-lo :E-hi E-hi)]
    (mapv (fn [{:keys [level energy]}]
            (let [e (exact level)]
              {:level level
               :numeric energy
               :exact e
               :rel-error (/ (Math/abs (- energy e)) (Math/abs e))}))
          levels)))

(defn spectrum
  "Solve a benchmark map for its bound-state energies (no closed form needed).
  Returns the vector of {:level :energy :nodes} maps."
  [{:keys [params search]} & {:keys [n-levels E-lo E-hi]}]
  (solve params
         :n-levels (or n-levels (:n-levels search) 4)
         :E-lo (or E-lo (:E-lo search) 0.0)
         :E-hi (or E-hi (:E-hi search) 200.0)))

(defn -main
  "Print benchmark validation for the closed-form potentials, plus the
  bound-state spectrum of the closed-form-less Woods-Saxon well."
  [& _]
  (println "Infinite square well (L=1, hbar=m=1):  E_n = n^2 pi^2 / 2")
  (doseq [row (benchmark-report (infinite-well {:L 1.0}))]
    (println (format "  level %d: numeric=%.6f exact=%.6f rel-err=%.2e"
                     (:level row) (:numeric row) (:exact row) (:rel-error row))))
  (println)
  (println "Harmonic oscillator (omega=1, hbar=m=1):  E_n = (n + 1/2)")
  (doseq [row (benchmark-report (harmonic {:omega 1.0}))]
    (println (format "  level %d: numeric=%.6f exact=%.6f rel-err=%.2e"
                     (:level row) (:numeric row) (:exact row) (:rel-error row))))
  (println)
  (println "Finite square well (V0=20, a=1):  transcendental closed form")
  (doseq [row (benchmark-report (finite-well {:V0 20.0 :a 1.0}))]
    (println (format "  level %d: numeric=%.6f exact=%.6f rel-err=%.2e"
                     (:level row) (:numeric row) (:exact row) (:rel-error row))))
  (println)
  (println "Woods-Saxon (V0=50, R=3, a=0.6):  NO closed form (grid-converged E_n)")
  (doseq [{:keys [level energy nodes]} (spectrum (woods-saxon {}))]
    (println (format "  level %d: E=%.6f  nodes=%d" level energy nodes))))
