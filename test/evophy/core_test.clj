(ns evophy.core-test
  (:require [clojure.test :refer :all]
            [evophy.core :refer :all]))

(deftest generate-data-smoke
  (testing "harmonic trajectory has expected length"
    (is (= 101 (count (vec (generate-data 1.0 1.0 1.0 0.0 0.1 100))))))

(deftest prediction-fitness-favors-true-harmonic-rates
  (let [data (vec (generate-data 1.0 1.0 1.0 0.0 0.1 20))
        ;; For m=k=1: dq/dt = p, dp/dt = -q
        true-rates {:dq 'p :dp '(- q)}
        flat-rates {:dq 0.5 :dp 0.5}
        true-fit (calculate-fitness (:dq true-rates) (:dp true-rates) data)
        flat-fit (calculate-fitness (:dq flat-rates) (:dp flat-rates) data)
        true-metrics (evaluate-predictions (:dq true-rates) (:dp true-rates) data)]
    (testing "true harmonic rates forecast better than constants"
      (is (> true-fit flat-fit)))
    (testing "true rates achieve small one-step error"
      (is (< (:mse true-metrics) 1e-4))))))
