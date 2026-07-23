(ns evophy.qm.matrix-test
  "Validate the matrix (finite-difference + Jacobi diagonalization) oracle
  against the same closed-form benchmarks as the shooting solver, and
  cross-check the two independent methods against each other. If they agree,
  neither is a fluke; if they disagree, at least one has a bug."
  (:require [clojure.test :refer [deftest is testing]]
            [evophy.qm.matrix :as m]
            [evophy.qm.schrodinger :as s]))

(defn- rel-err [a b] (/ (Math/abs (- a b)) (Math/abs b)))

;; O(n^3) diagonalization -- keep the grid modest for test speed.
(def ^:private N 201)

(deftest infinite-well-spectrum
  (testing "E_n = n^2 pi^2 / (2 m L^2) for L=1, hbar=m=1"
    (let [{:keys [params exact]} (s/infinite-well {:L 1.0 :n N})
          levels (m/solve params :n-levels 4)]
      (is (= 4 (count levels)) "found four bound levels")
      (doseq [{:keys [level energy]} levels]
        (is (< (rel-err energy (exact level)) 1e-2)
            (format "level %d numeric %.5f vs exact %.5f" level energy (exact level)))))))

(deftest harmonic-spectrum
  (testing "E_n = (n + 1/2) for omega=1, hbar=m=1"
    (let [{:keys [params exact]} (s/harmonic {:omega 1.0 :x-max 12.0 :n N})
          levels (m/solve params :n-levels 5)]
      (is (= 5 (count levels)) "found five bound levels")
      (doseq [{:keys [level energy]} levels]
        (is (< (rel-err energy (exact level)) 1e-2)
            (format "level %d numeric %.5f vs exact %.5f" level energy (exact level)))))))

(deftest node-count-matches-level
  (testing "the k-th bound state has exactly k interior nodes (discrete oscillation theorem)"
    (let [{:keys [params]} (s/harmonic {:omega 1.0 :x-max 12.0 :n N})
          levels (m/solve params :n-levels 4)]
      (doseq [{:keys [level nodes]} levels]
        (is (= level nodes)
            (format "level %d should have %d nodes, got %d" level level nodes))))))

(deftest finite-well-spectrum
  (testing "finite square well energies match the transcendental closed form"
    ;; The finite-difference stencil is O(h^2) (vs Numerov's O(h^4)) and, like
    ;; the shooting solver, degrades further at the well's discontinuous edge
    ;; -- tolerance is looser here than in schrodinger_test's equivalent case
    ;; for that reason, not because fewer levels converge.
    (let [bench (s/finite-well {:V0 20.0 :a 1.0 :x-max 8.0 :n N})
          rows  (m/benchmark-report bench)
          exact-levels (s/finite-well-levels 20.0 1.0)]
      (is (= (count exact-levels) (count rows))
          "solver finds exactly the transcendental bound-state count")
      (doseq [{:keys [level numeric exact rel-error]} rows
              :when (> (Math/abs exact) 0.5)]
        (is (< rel-error 5e-2)
            (format "level %d numeric %.5f vs exact %.5f" level numeric exact))))))

(deftest wavefunction-normalized
  (testing "returned wavefunctions are unit-normalized on the grid"
    (let [{:keys [params]} (s/harmonic {:omega 1.0 :x-max 12.0 :n N})
          h (s/step-size params)
          {:keys [psi]} (first (m/solve params :n-levels 1))
          norm2 (reduce + (map (fn [v] (* v v h)) psi))]
      (is (< (Math/abs (- norm2 1.0)) 1e-6)
          (format "norm^2 = %.5f" norm2)))))

(deftest agrees-with-shooting-solver
  (testing "matrix and Numerov+shooting agree on Woods-Saxon (no closed form to check either against)"
    (let [bench (s/woods-saxon {:n N})
          mat   (m/spectrum bench :n-levels 5)
          shoot (s/spectrum bench :n-levels 5)]
      (is (= (count mat) (count shoot)) "both methods find the same number of levels")
      (doseq [[a b] (map vector mat shoot)]
        (is (< (Math/abs (- (:energy a) (:energy b))) 1e-1)
            (format "level %d: matrix %.6f vs shooting %.6f" (:level a) (:energy a) (:energy b)))))))
