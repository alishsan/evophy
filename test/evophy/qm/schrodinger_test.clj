(ns evophy.qm.schrodinger-test
  "Validate the Numerov shooting solver against the two closed-form benchmarks.
  These are the QM analogue of checking the symplectic integrator against
  known Kepler orbits: if the ground-truth oracle is wrong, nothing built on
  top can be trusted."
  (:require [clojure.test :refer [deftest is testing]]
            [evophy.qm.schrodinger :as s]))

(defn- rel-err [a b] (/ (Math/abs (- a b)) (Math/abs b)))

(deftest infinite-well-spectrum
  (testing "E_n = n^2 pi^2 / (2 m L^2) for L=1, hbar=m=1"
    (let [{:keys [params exact]} (s/infinite-well {:L 1.0 :n 4001})
          levels (s/solve params :n-levels 4 :E-hi 250.0)]
      (is (= 4 (count levels)) "found four bound levels")
      (doseq [{:keys [level energy]} levels]
        (is (< (rel-err energy (exact level)) 1e-3)
            (format "level %d numeric %.5f vs exact %.5f" level energy (exact level)))))))

(deftest harmonic-spectrum
  (testing "E_n = (n + 1/2) for omega=1, hbar=m=1"
    (let [{:keys [params exact]} (s/harmonic {:omega 1.0 :x-max 12.0 :n 4001})
          levels (s/solve params :n-levels 5 :E-hi 20.0)]
      (is (= 5 (count levels)) "found five bound levels")
      (doseq [{:keys [level energy]} levels]
        (is (< (rel-err energy (exact level)) 1e-3)
            (format "level %d numeric %.5f vs exact %.5f" level energy (exact level)))))))

(deftest node-count-matches-level
  (testing "the k-th bound state has exactly k interior nodes"
    (let [{:keys [params]} (s/harmonic {:omega 1.0 :x-max 12.0 :n 4001})
          levels (s/solve params :n-levels 4 :E-hi 20.0)]
      (doseq [{:keys [level nodes]} levels]
        (is (= level nodes)
            (format "level %d should have %d nodes, got %d" level level nodes))))))

(deftest finite-well-spectrum
  (testing "finite square well energies match the transcendental closed form"
    ;; x-max=8, n=4001 puts the well edges (|x|=1) exactly on grid points.
    ;; Tolerance is 2e-2, not machine-tight: the well's discontinuous potential
    ;; step degrades Numerov from O(h^4) to O(h) at the edge, so this benchmark
    ;; is genuinely grid-limited (unlike the smooth harmonic).
    (let [bench (s/finite-well {:V0 20.0 :a 1.0 :x-max 8.0 :n 4001})
          rows  (s/benchmark-report bench)
          exact-levels (s/finite-well-levels 20.0 1.0)]
      (is (= (count exact-levels) (count rows))
          "solver finds exactly the transcendental bound-state count")
      ;; Deep, well-resolved states must agree. The near-threshold state
      ;; (|E| < 0.5) has an essentially unbounded tail a finite box can't
      ;; represent, so it is intentionally exempt from the accuracy check.
      (doseq [{:keys [level numeric exact rel-error]} rows
              :when (> (Math/abs exact) 0.5)]
        (is (< rel-error 2e-2)
            (format "level %d numeric %.5f vs exact %.5f" level numeric exact))))))

(deftest woods-saxon-grid-convergence
  (testing "closed-form-less Woods-Saxon well: consistent, ordered, grid-converged"
    (let [coarse (s/spectrum (s/woods-saxon {:n 1501}) :n-levels 6)
          fine   (s/spectrum (s/woods-saxon {:n 3001}) :n-levels 6)]
      (is (= 6 (count coarse)) "coarse grid finds six bound states")
      (is (= 6 (count fine))   "fine grid finds six bound states")
      (doseq [{:keys [level nodes]} fine]
        (is (= level nodes)
            (format "level %d should have %d nodes, got %d" level level nodes)))
      (is (apply < (map :energy fine)) "energies strictly increasing")
      (doseq [[c f] (map vector coarse fine)]
        (is (< (Math/abs (- (:energy c) (:energy f))) 5e-2)
            (format "level %d: coarse %.4f vs fine %.4f grid" (:level f) (:energy c) (:energy f)))))))

(deftest wavefunction-normalized
  (testing "returned wavefunctions are unit-normalized on the grid"
    (let [{:keys [params]} (s/harmonic {:omega 1.0 :x-max 12.0 :n 4001})
          h (s/step-size params)
          {:keys [psi]} (first (s/solve params :n-levels 1 :E-hi 20.0))
          norm2 (reduce + (map (fn [v] (* v v h)) psi))]
      (is (< (Math/abs (- norm2 1.0)) 1e-2)
          (format "norm^2 = %.5f" norm2)))))
