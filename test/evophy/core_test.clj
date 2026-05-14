(ns evophy.core-test
  (:require [clojure.test :refer :all]
            [evophy.core :refer :all]))

(deftest generate-data-smoke
  (testing "harmonic trajectory has expected length"
    (is (= 101 (count (vec (generate-data 1.0 1.0 1.0 0.0 0.1 100)))))))

(deftest calculate-fitness-distinguishes-literals-from-state-laws
  (let [data (vec (generate-data 1.0 1.0 1.0 0.0 0.1 20))
        h-expr '(+ (* 0.5 (e/square p)) (* 0.5 (e/square q)))
        h-fit (calculate-fitness h-expr data)]
    (testing "numeric literals are flat in (q,p), zero fitness"
      (is (zero? (calculate-fitness 1.0 data)))
      (is (zero? (calculate-fitness 0.5 data))))
    (testing "Hamiltonian-shaped expression scores and beats constants"
      (is (pos? h-fit))
      (is (> h-fit (calculate-fitness 0.5 data))))))
