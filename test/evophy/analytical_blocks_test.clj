(ns evophy.analytical-blocks-test
  (:require [clojure.test :refer :all]
            [evophy.analytical-blocks :as blocks]
            [evophy.core :as core :refer :all]
            [evophy.symbolic :as sym]))

(deftest catalog-loads-with-emmy-validation
  (is (= 6 (blocks/catalog-block-count)))
  (is (some #(= :circle-omega (:id %)) blocks/validated-catalog))
  (is (some #(= :ellipse-omega-L (:id %)) blocks/validated-catalog))
  (is (some #(= :ellipse-conic (:id %)) blocks/validated-catalog))
  (is (some #(= :hyperbola-conic (:id %)) blocks/validated-catalog))
  (is (not (some #(= :taylor-bound (:id %)) blocks/validated-catalog)))
  (is (not (some #(= :hyperbola (:id %)) blocks/validated-catalog))))

(deftest kepler-proxy-satisfies-expr-uses-t
  (is (expr-uses-t? 'kepler-u))
  (is (expr-uses-t? '(e/cos kepler-u)))
  (is (expr-uses-t? '(e/sinh kepler-F)))
  (is (not (expr-uses-t? '(+ q0x p0x)))))

(deftest conic-catalog-de-fitness-beats-axis-mix
  (let [ds (scenarios->datasets default-scenarios)]
    (binding [*both-regimes?* true *template-unbound-arms?* true]
      (let [conic-fit (#'core/calculate-analytical-de-driven-fitness
                       (#'core/wrap-catalog-entry-strict (blocks/catalog-entry :ellipse-conic))
                       ds)
            axis-fit (#'core/calculate-analytical-de-driven-fitness
                      (#'core/wrap-catalog-entry-strict (blocks/catalog-entry :ellipse-axis-mix))
                      ds)]
        (is (pos? conic-fit))
        (is (> conic-fit axis-fit))))))

(deftest kepler-conic-full-regime-fitness
  (let [ds (scenarios->datasets default-scenarios)]
    (binding [*both-regimes?* true
              *template-unbound-arms?* true
              *template-conic-unbound?* true
              *analytical-blocks?* true]
      (let [kepler (#'core/wrap-kepler-conic-strict (blocks/kepler-conic-entry))
            ell-only (#'core/wrap-catalog-entry-strict (blocks/catalog-entry :ellipse-conic))
            k-fit (#'core/calculate-analytical-de-driven-fitness kepler ds)
            e-fit (#'core/calculate-analytical-de-driven-fitness ell-only ds)]
        (is (analytical-strict-energy-branches? kepler))
        (is (#'core/expr-uses-op? (unbound-arm-expr :qx-expr) 'e/sinh))
        (is (> k-fit e-fit))
        (is (> k-fit 0.95))))))

(deftest conic-catalog-matches-kepler
  (let [ell (blocks/catalog-entry :ellipse-conic)
        hyp (blocks/catalog-entry :hyperbola-conic)
        compile (fn [entry]
                  (#'core/compile-analytical-fns
                   (into {:strategy :analytical} (:slots entry))
                   :cartesian))]
    (is (some? ell))
    (is (some? hyp))
    (let [fns (compile ell)
          m 1.0 a 1.0 q0x 1.2 q0y 0.0 p0x 0.0 p0y 0.9 t 0.5
          ks (#'core/kepler-state-at-t m a q0x q0y p0x p0y t)]
      (is (some? ks))
      (is (< (Math/abs (- ((:qx fns) t q0x q0y p0x p0y m a) (:qx ks))) 1e-9))
      (is (< (Math/abs (- ((:qy fns) t q0x q0y p0x p0y m a) (:qy ks))) 1e-9)))
    (let [fns (compile hyp)
          m 1.0 a 1.0 q0x 2.5 q0y 0.3 p0x -0.4 p0y 0.8 t 0.5
          ks (#'core/kepler-state-at-t m a q0x q0y p0x p0y t)]
      (is (some? ks))
      (is (< (Math/abs (- ((:qx fns) t q0x q0y p0x p0y m a) (:qx ks))) 1e-9))
      (is (< (Math/abs (- ((:py fns) t q0x q0y p0x p0y m a) (:py ks))) 1e-9)))))

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
