(ns evophy.analytical-blocks-test
  (:require [clojure.test :refer :all]
            [evophy.analytical-blocks :as blocks]
            [evophy.core :refer :all]
            [evophy.symbolic :as sym]))

(deftest catalog-loads-with-emmy-validation
  (is (>= (blocks/catalog-block-count) 5))
  (is (some #(= :circle-omega (:id %)) blocks/validated-catalog))
  (is (some #(= :ellipse-omega-L (:id %)) blocks/validated-catalog)))

(deftest catalog-circle-block-compiles
  (let [law (blocks/random-catalog-law)]
    (binding [*both-regimes?* true
              *template-unbound-arms?* true]
      (let [wrapped (template-unbound-arms-law
                     (into law
                           (map (fn [k] [k (slot-with-bound-arm k (get law k))])
                                [:qx-expr :qy-expr :px-expr :py-expr])))]
        (is (analytical-strict-energy-branches? wrapped))
        (is (genome-valid? wrapped))))))

(deftest symbolic-validates-taylor
  (is (sym/symbolic-validate-expr?
       '(+ q0x (* (e/div p0x m) t)))))

(deftest graft-catalog-pair-preserves-strict-branches
  (binding [*both-regimes?* true
            *template-unbound-arms?* true
            *analytical-blocks?* true]
    (let [branch (fn [ex] (list 'e/if '(neg? energy) ex ex))
          taylor {:strategy :analytical :domain :any
                  :qx-expr '(+ q0x (* (e/div p0x m) t))
                  :qy-expr '(+ q0y (* (e/div p0y m) t))
                  :px-expr '(+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
                  :py-expr '(+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))}
          parent (into taylor
                       (map (fn [[k ex]] [k (branch ex)])
                            (select-keys taylor [:qx-expr :qy-expr :px-expr :py-expr])))
          entry (blocks/catalog-entry :circle-omega)
          grafted (blocks/graft-bound-pair parent entry [:qx-expr :px-expr])
          wrapped (template-unbound-arms-law
                   (into grafted
                         (map (fn [k] [k (slot-with-bound-arm k (get grafted k))])
                              [:qx-expr :px-expr])))]
      (is (analytical-strict-energy-branches? wrapped))
      (is (not= (get parent :qx-expr) (get wrapped :qx-expr))))))
