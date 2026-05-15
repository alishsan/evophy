(ns evophy.core-test
  (:require [clojure.java.io :as io]
            [clojure.test :refer :all]
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

(deftest multi-scenario-fitness-generalizes-true-rates
  (let [datasets (scenarios->datasets (take 3 default-scenarios))
        true-rates {:dq 'p :dp '(- q)}
        flat-rates {:dq 0.5 :dp 0.5}
        true-fit (calculate-fitness-scenarios (:dq true-rates) (:dp true-rates) datasets)
        flat-fit (calculate-fitness-scenarios (:dq flat-rates) (:dp flat-rates) datasets)]
    (testing "true harmonic rates beat constants averaged over scenarios"
      (is (> true-fit flat-fit)))
    (testing "each scenario dataset is non-empty"
      (is (every? #(pos? (count %)) datasets)))))

(deftest population-persistence-roundtrip
  (let [path (str (io/file (System/getProperty "java.io.tmpdir")
                           (str "evophy-pop-" (System/nanoTime) ".edn")))
        pop [{:dq 'p :dp '(- q)} {:dq 'q :dp 'p}]]
    (try
      (save-population! path pop :generations-run 12 :population-size 2)
      (is (= 12 (:generations-run (load-population path))))
      (is (= pop (:population (load-population path))))
      (finally
        (.delete (io/file path)))))

(deftest resolve-initial-population-respects-fresh
  (let [path (str (io/file (System/getProperty "java.io.tmpdir")
                           (str "evophy-pop-" (System/nanoTime) ".edn")))
        saved [{:dq 1 :dp 2}]]
    (try
      (save-population! path saved :generations-run 5 :population-size 1)
      (is (:resumed? (resolve-initial-population {:fresh? false
                                                   :path path
                                                   :population-size 1})))
      (is (not (:resumed? (resolve-initial-population {:fresh? true
                                                       :path path
                                                       :population-size 1}))))
      (finally
        (.delete (io/file path)))))))
