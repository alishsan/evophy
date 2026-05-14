(ns evophy.core-test
  (:require [clojure.test :refer :all]
            [evophy.core :refer :all]))

(deftest generate-data-smoke
  (testing "harmonic trajectory has expected length"
    (is (= 101 (count (vec (generate-data 1.0 1.0 1.0 0.0 0.1 100)))))))
