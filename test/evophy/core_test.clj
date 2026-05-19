(ns evophy.core-test
  (:require [clojure.java.io :as io]
            [clojure.test :refer :all]
            [evophy.core :refer :all]))

(deftest generate-data-smoke
  (testing "2D gravitational trajectory has expected length"
    (is (= 101 (count (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 100)))))))

(deftest gravitational-trajectory-varies-in-time
  (let [data (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 40))
        qs (map :qx data)]
    (is (pos? (count data)))
    (is (not= (first qs) (last qs)))))

(deftest calculate-fitness-smoke-gravitational
  (let [data (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 20))
        ind (random-analytical-individual)]
    (is (number? (calculate-fitness ind data)))))

(deftest t-only-analytical-laws-score-zero
  (let [data (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 20))
        cheat {:strategy :analytical
               :qx-expr '(e/sin t) :qy-expr '(e/sin t)
               :px-expr '(e/cos t) :py-expr '(e/cos t)}]
    (testing "laws that ignore ICs are rejected"
      (is (zero? (calculate-fitness cheat data))))))

(deftest analytical-missing-ic-symbol-scores-zero
  (let [data (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 20))
        cheat {:strategy :analytical
               :qx-expr '(e/sin t)
               :qy-expr '(e/sin t)
               :px-expr '(e/sin t)
               :py-expr '(e/sin t)}]
    (testing "laws that never use p0x are rejected"
      (is (zero? (calculate-fitness cheat data))))))

(deftest analytical-exprs-must-use-time
  (let [data (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 20))
        static {:strategy :analytical
                :qx-expr '(e/square (* p0x q0x))
                :qy-expr '(e/square (* p0y q0y))
                :px-expr '(* t p0x q0x)
                :py-expr '(* t p0y q0y)}]
    (testing "qx with no t is rejected"
      (is (zero? (calculate-fitness static data))))))

(def ^:private r0cubed
  "r₀³ with q0x/q0y only (no syntax-quote — avoids namespaced symbols in compiled fns)."
  '(* (e/sqrt (+ (e/square q0x) (e/square q0y)))
      (e/square (e/sqrt (+ (e/square q0x) (e/square q0y))))))

(deftest differential-fitness-favors-true-gravitational-rates
  "At t0 with m=α=1: q̇x=p0x, q̇y=p0y, ṗx=−q0x/r₀³, ṗy=−q0y/r₀³."
  (let [data (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 20))
        true-ind {:strategy :differential
                  :dqx-expr 'p0x
                  :dqy-expr 'p0y
                  :dpx-expr (list '* -1.0 (list 'e/div 'q0x r0cubed))
                  :dpy-expr (list '* -1.0 (list 'e/div 'q0y r0cubed))}
        flat-ind {:strategy :differential
                  :dqx-expr 'q0x :dqy-expr 'q0y :dpx-expr 'p0x :dpy-expr 'p0y}]
    (testing "true 2D gravitational rates beat wrong linear rates"
      (is (pos? (calculate-fitness true-ind data)))
      (is (> (calculate-fitness true-ind data)
             (calculate-fitness flat-ind data))))))

(deftest multi-scenario-fitness-generalizes-true-rates
  (let [datasets (scenarios->datasets (take 3 default-scenarios))
        true-ind {:strategy :differential
                  :dqx-expr 'p0x
                  :dqy-expr 'p0y
                  :dpx-expr '(* -1.0 (e/div q0x r0cubed))
                  :dpy-expr '(* -1.0 (e/div q0y r0cubed))}
        flat-ind {:strategy :differential
                  :dqx-expr 'q0x :dqy-expr 'q0y :dpx-expr 'p0x :dpy-expr 'p0y}
        true-fit (calculate-fitness-scenarios true-ind datasets)
        flat-fit (calculate-fitness-scenarios flat-ind datasets)]
    (testing "true rates beat constants on worst scenario"
      (is (> true-fit flat-fit)))
    (testing "each scenario dataset is non-empty"
      (is (every? #(pos? (count %)) datasets)))))

(deftest population-persistence-roundtrip
  (let [path (str (io/file (System/getProperty "java.io.tmpdir")
                           (str "evophy-pop-" (System/nanoTime) ".edn")))
        pop [{:strategy :analytical
              :qx-expr 'p0x :qy-expr 'p0y :px-expr 't :py-expr 't}
             {:strategy :differential
              :dqx-expr 'p0x :dqy-expr 'p0y
              :dpx-expr '(* -1.0 (e/div q0x (e/square q0x)))
              :dpy-expr '(* -1.0 (e/div q0y (e/square q0y)))}]]
    (try
      (save-population! path pop :generations-run 12 :population-size 2)
      (is (= 12 (:generations-run (load-population path))))
      (is (= pop (:population (load-population path))))
      (finally
        (.delete (io/file path))))))

(deftest resolve-initial-population-respects-fresh
  (let [path (str (io/file (System/getProperty "java.io.tmpdir")
                           (str "evophy-pop-" (System/nanoTime) ".edn")))
        saved [{:strategy :analytical
                :qx-expr 1 :qy-expr 2 :px-expr 3 :py-expr 4}]]
    (try
      (save-population! path saved :generations-run 5 :population-size 1)
      (is (:resumed? (resolve-initial-population {:fresh? false
                                                   :path path
                                                   :population-size 1})))
      (is (not (:resumed? (resolve-initial-population {:fresh? true
                                                       :path path
                                                       :population-size 1}))))
      (finally
        (.delete (io/file path))))))
