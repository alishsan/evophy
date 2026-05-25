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
