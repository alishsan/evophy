(ns evophy.mcts-test
  (:require [clojure.test :refer :all]
            [evophy.core :refer :all]
            [evophy.mcts :as mcts]))

(deftest mcts-produces-valid-analytical-genome
  (let [datasets (scenarios->datasets (take 2 default-scenarios))
        ind (mcts/search-analytical-individual datasets 64)]
    (is (= :analytical (:strategy ind)))
    (is (genome-valid? ind))
    (is (number? (calculate-fitness-scenarios ind datasets)))))
