(ns evophy.mcts-test
  (:require [clojure.test :refer :all]
            [evophy.core :refer :all]
            [evophy.mcts :as mcts]))

(deftest mcts-produces-valid-analytical-genome
  (let [datasets (scenarios->datasets (take 2 default-scenarios))
        _ (mcts/clear-stop!)
        ind (mcts/search-analytical-individual datasets 64)]
    (is (= :analytical (:strategy ind)))
    (is (genome-valid? ind))
    (is (number? (calculate-fitness-scenarios ind datasets)))))

(deftest mcts-stops-early-with-custom-should-stop
  (let [datasets (scenarios->datasets (take 1 default-scenarios))
        _ (mcts/clear-stop!)
        n (atom 0)
        ind (mcts/search-analytical-individual datasets 1000
               :should-stop? #(>= (swap! n inc) 4))
        stats (mcts/last-search-result)]
    (is (= :analytical (:strategy ind)))
    (is (:stopped? stats))
    (is (<= (:simulations-run stats) 4))
    (is (pos? (:simulations-run stats)))))

(deftest adversarial-repair-target-finds-weak-slot
  (let [circle (scenario-data (first default-scenarios))
        hyper  (scenario-data (last default-scenarios))
        taylor {:strategy :analytical :domain :unbound
                :qx-expr '(+ q0x (* (e/div p0x m) t))
                :qy-expr '(+ q0y (* (e/div p0y m) t))
                :px-expr '(+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
                :py-expr '(+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))}
        datasets (scenarios->datasets default-scenarios)
        target (adversarial-repair-target taylor datasets)]
    (is (map? target))
    (is (contains? target :dataset))
    (is (contains? #{:qx-expr :qy-expr :px-expr :py-expr} (:expr-key target)))))

(deftest mcts-repair-produces-valid-analytical-genome
  (let [datasets (scenarios->datasets (take 3 default-scenarios))
        branch (fn [ex] (list 'e/if '(neg? energy) ex ex))
        taylor {:strategy :analytical :domain :any
                :qx-expr '(+ q0x (* (e/div p0x m) t))
                :qy-expr '(+ q0y (* (e/div p0y m) t))
                :px-expr '(+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
                :py-expr '(+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))}
        seed (into taylor
                   (map (fn [[k ex]] [k (branch ex)])
                        (select-keys taylor [:qx-expr :qy-expr :px-expr :py-expr])))
        _ (mcts/clear-stop!)
        repaired (binding [*both-regimes?* true]
                   (mcts/search-repair-individual seed datasets 32))]
    (is (genome-valid? repaired))
    (is (analytical-strict-energy-branches? repaired))
    (is (number? (calculate-fitness-scenarios repaired datasets :evaluation :de-driven)))
    (let [stats (mcts/last-search-result)]
      (is (:repair? stats))
      (is (map? (:repair-target stats))))))

(deftest mcts-mutate-wiring-produces-valid-child
  (let [datasets (scenarios->datasets (take 2 default-scenarios))
        branch (fn [ex] (list 'e/if '(neg? energy) ex ex))
        taylor {:strategy :analytical :domain :any
                :qx-expr '(+ q0x (* (e/div p0x m) t))
                :qy-expr '(+ q0y (* (e/div p0y m) t))
                :px-expr '(+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
                :py-expr '(+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))}
        parent (into taylor
                       (map (fn [[k ex]] [k (branch ex)])
                            (select-keys taylor [:qx-expr :qy-expr :px-expr :py-expr])))
        _ (mcts/clear-stop!)
        child (binding [*mcts-mutate?* true
                        *mcts-mutate-rate* 1.0
                        *mcts-mutate-datasets* datasets
                        *mcts-mutate-simulations* 12
                        *both-regimes?* true
                        *guess-mutations?* false
                        *stagnation-escape?* false]
                (with-redefs [rand (constantly 0.99)]
                  (mutate-individual parent)))]
    (is (genome-valid? child))
    (is (analytical-strict-energy-branches? child))
    (is (= :analytical (:strategy child)))))
