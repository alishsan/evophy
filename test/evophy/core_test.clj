(ns evophy.core-test
  (:require [clojure.java.io :as io]
            [clojure.test :refer :all]
            [evophy.core :refer :all]))

(defn- dataset
  "Wrap raw trajectory with scenario parameters for fitness APIs."
  [m alpha data]
  {:m m :alpha alpha :data data})

(def ^:private r-cubed
  "r³ from current position (qx, qy)."
  '(* (e/sqrt (+ (e/square qx) (e/square qy)))
      (e/square (e/sqrt (+ (e/square qx) (e/square qy))))))

(def ^:private true-grav-rates
  "Hamiltonian rates: q̇ = p/m, ṗ = −α q/r³."
  {:strategy :differential
   :dqx-expr '(e/div px m)
   :dqy-expr '(e/div py m)
   :dpx-expr (list '* -1.0 (list '* 'alpha (list 'e/div 'qx r-cubed)))
   :dpy-expr (list '* -1.0 (list '* 'alpha (list 'e/div 'qy r-cubed)))})

(def ^:private wrong-grav-rates
  "Uses m and α but wrong dependence on state."
  {:strategy :differential
   :dqx-expr '(* qx m)
   :dqy-expr '(* qy m)
   :dpx-expr '(* px alpha)
   :dpy-expr '(* py alpha)})

(deftest generate-data-smoke
  (testing "2D gravitational trajectory has expected length"
    (is (= 101 (count (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 100)))))))

(deftest gravitational-trajectory-varies-in-time
  (let [data (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 40))
        qs (map :qx data)]
    (is (pos? (count data)))
    (is (not= (first qs) (last qs)))))

(deftest calculate-fitness-smoke-gravitational
  (let [ds (dataset 1.0 1.0 (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 20)))
        ind (random-analytical-individual)]
    (is (number? (calculate-fitness ind ds)))))

(deftest t-only-analytical-laws-score-zero
  (let [ds (dataset 1.0 1.0 (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 20)))
        cheat {:strategy :analytical
               :qx-expr '(e/sin t) :qy-expr '(e/sin t)
               :px-expr '(e/cos t) :py-expr '(e/cos t)}]
    (testing "laws that ignore ICs are rejected"
      (is (zero? (calculate-fitness cheat ds))))))

(deftest analytical-missing-ic-symbol-scores-zero
  (let [ds (dataset 1.0 1.0 (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 20)))
        cheat {:strategy :analytical
               :qx-expr '(e/sin t)
               :qy-expr '(e/sin t)
               :px-expr '(e/sin t)
               :py-expr '(e/sin t)}]
    (testing "laws that never use p0x are rejected"
      (is (zero? (calculate-fitness cheat ds))))))

(deftest analytical-exprs-must-use-time
  (let [ds (dataset 1.0 1.0 (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 20)))
        static {:strategy :analytical
                :qx-expr '(e/square (* p0x q0x))
                :qy-expr '(e/square (* p0y q0y))
                :px-expr '(* t p0x q0x m)
                :py-expr '(* t p0y q0y alpha)}]
    (testing "qx with no t is rejected"
      (is (zero? (calculate-fitness static ds))))))

(deftest differential-fitness-favors-true-gravitational-rates
  (let [ds (dataset 1.0 1.0 (vec (generate-data 1.0 1.0 2.5 0.0 0.0 0.35 0.05 20)))]
    (testing "true parameterized rates beat wrong rates"
      (is (pos? (calculate-fitness true-grav-rates ds)))
      (is (> (calculate-fitness true-grav-rates ds)
             (calculate-fitness wrong-grav-rates ds))))))

(deftest multi-scenario-fitness-generalizes-true-rates
  (let [datasets (scenarios->datasets default-scenarios)
        true-fit (calculate-fitness-scenarios true-grav-rates datasets)
        wrong-fit (calculate-fitness-scenarios wrong-grav-rates datasets)]
    (testing "true rates beat wrong rates on worst scenario"
      (is (> true-fit wrong-fit)))
    (testing "each scenario dataset is non-empty"
      (is (every? #(pos? (count (:data %))) datasets)))))

(deftest random-scenario-respects-radius-floor
  (let [{:keys [q0x q0y]} (first (sample-random-scenarios 1 :seed 42))]
    (is (>= (Math/sqrt (+ (* q0x q0x) (* q0y q0y))) 1.0))))

(deftest aggregate-percentile-between-min-and-mean
  (let [fits [0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]]
    (is (>= (aggregate-scenario-fitness fits :aggregate :percentile :percentile 0.2)
            (aggregate-scenario-fitness fits :aggregate :min)))
    (is (<= (aggregate-scenario-fitness fits :aggregate :percentile :percentile 0.2)
            (/ (reduce + fits) (count fits))))))

(deftest true-rates-beat-wrong-on-random-scenario-batch
  (let [ctx (make-fitness-context :mode :random :sample-count 6 :seed 99 :aggregate :min)
        datasets (datasets-for-fitness-context ctx)]
    (is (= 6 (count datasets)))
    (is (> (calculate-fitness-scenarios true-grav-rates datasets)
           (calculate-fitness-scenarios wrong-grav-rates datasets)))))

(deftest datasets-for-context-new-batch-each-generation
  (let [ctx (make-fitness-context :mode :random :sample-count 8 :seed 1)
        a (datasets-for-fitness-context ctx :generation 0)
        b (datasets-for-fitness-context ctx :generation 1)]
    (is (not= (map :q0x a) (map :q0x b)))))

(deftest true-rates-fit-heavy-m-scenario
  (let [ds (scenario-data {:m 2.0 :alpha 1.0 :q0x 2.2 :q0y 1.0 :p0x 0.1 :p0y 0.25 :dt 0.04 :steps 40})]
    (is (pos? (calculate-fitness true-grav-rates ds)))))

(deftest simplifier-eliminates-dead-algebra
  (is (= 0 (#'evophy.core/simplify-expr '(* (- qx qx) alpha))))
  (is (= 'qx (#'evophy.core/simplify-expr '(+ qx 0))))
  (is (= 'px (#'evophy.core/simplify-expr '(* 1 px)))))

(deftest behavior-key-collapses-equivalent-expressions
  (let [probes (build-behavior-probes (scenarios->datasets (take 2 default-scenarios)))
        variant (assoc true-grav-rates :dpx-expr (list '+ (list '* -1.0 (list '* 'alpha (list 'e/div 'qx r-cubed))) (list '- 'qx 'qx)))]
    (is (= (individual-behavior-key true-grav-rates probes)
           (individual-behavior-key variant probes)))))

(deftest distinct-by-behavior-keeps-fewer-when-synonyms-would-dup-genome
  (let [probes (build-behavior-probes (scenarios->datasets (take 2 default-scenarios)))
        variant (assoc true-grav-rates :dpx-expr (list '+ (list '* -1.0 (list '* 'alpha (list 'e/div 'qx r-cubed))) (list '- 'qx 'qx)))
        ranked [(assoc variant :fitness 99.0) (assoc true-grav-rates :fitness 100.0)]]
    (is (= 1 (count (take-distinct-by-behavior 5 ranked probes))))))

(deftest mutate-produces-printable-expressions
  (let [tree '(+ (e/square qx) (* py m))
        mutated (mutate tree differential-vars)]
    (is (list? mutated))
    (is (not (instance? clojure.lang.LazySeq mutated)))))

(deftest individual->equations-renders-laws
  (let [ind (assoc true-grav-rates
              :dqx-expr 'p0x :dqy-expr 'p0y
              :dpx-expr '(* -1.0 q0x) :dpy-expr '(* -1.0 q0y))
        plain (individual->equations ind :scenario-id :x-orbit)
        latex (individual->equations ind :scenario-id :x-orbit :format :latex)]
    (is (= :differential (:strategy plain)))
    (is (clojure.string/includes? (:dqx (:equations plain)) "dq_x/dt"))
    (is (clojure.string/includes? (:dqx (:equations latex)) "\\dot{q}_x"))
    (is (map? (:metrics plain)))
    (is (contains? (:sample plain) :rates))))

(deftest population->motif-report-finds-shared-subtrees
  (let [pop [{:strategy :analytical
              :qx-expr '(+ (e/cos t) q0y)
              :qy-expr '(e/cos (e/square t))
              :px-expr 'p0x
              :py-expr 't}
             {:strategy :analytical
              :qx-expr '(+ (e/cos t) q0y)
              :qy-expr '(e/sin t)
              :px-expr 'p0x
              :py-expr 't}]
        report (population->motif-report pop :elite-frac 1.0 :min-share 0.5 :top 10)]
    (is (= 2 (:n-population report)))
    (is (pos? (count (:motifs report))))
    (is (some #(= (:signature %) "(+ (e/cos t) q0y)")
              (:motifs report)))))

(deftest individual->equations-analytical-sample
  (let [ind {:strategy :analytical
             :qx-expr '(+ q0x (* p0x t))
             :qy-expr '(+ q0y (* p0y t))
             :px-expr 'p0x
             :py-expr 'p0y}
        eq (individual->equations ind :scenario-id :x-orbit
                                  :sample-t 0.05)]
    (is (clojure.string/includes? (:qx (:equations eq)) "q_x(t)"))
    (is (= 0.05 (:t (:sample eq))))
    (is (number? (get-in eq [:sample :predicted :qx])))))

(deftest population-persistence-roundtrip
  (let [path (str (io/file (System/getProperty "java.io.tmpdir")
                           (str "evophy-pop-" (System/nanoTime) ".edn")))
        pop [{:strategy :analytical
              :qx-expr 'p0x :qy-expr 'p0y :px-expr 't :py-expr 't}
             true-grav-rates]]
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
