(ns evophy.coords-test
  (:require [clojure.test :refer :all]
            [evophy.coords :as coords]
            [evophy.core :as core]))

(deftest circle-scenario-prefers-polar
  (let [circle (first (filter #(= (:id %) :circle) core/default-scenarios))
        ds     (core/scenario-data circle)]
    (is (= :polar (:chart ds)))
    (is (every? #(contains? % :r) (:data ds)))))

(deftest tilted-ellipse-stays-cartesian
  (let [sc (first (filter #(= (:id %) :ellipse-tilted) core/default-scenarios))
        ds (core/scenario-data sc)]
    (is (= :cartesian (:chart ds)))))

(deftest polar-angular-momentum-is-ptheta
  (let [{:keys [ptheta]} (coords/cartesian->polar {:qx 2.0 :qy 0.0 :px 0.0 :py 1.414})]
    (is (< (Math/abs (- ptheta 2.828)) 0.01))))

(deftest effective-chart-defaults-from-scenario
  (let [circle (first (filter #(= (:id %) :circle) core/default-scenarios))
        ds     (core/scenario-data circle)
        laws   (:laws (first core/de-driven-composite-seeds))
        analytical (assoc (first laws) :strategy :analytical)
        conserved  (assoc (second laws) :strategy :conserved)]
    (is (= :polar (coords/effective-chart analytical (:chart ds))))
    (is (= :polar (coords/effective-chart conserved (:chart ds))))
    (is (pos? (core/calculate-conserved-fitness conserved ds)))))

(deftest cartesian-law-scores-zero-on-polar-scenario
  (let [ds  (core/scenario-data (first (filter #(= (:id %) :circle) core/default-scenarios)))
        law {:strategy :analytical
             :qx-expr '(+ q0x (* (e/div p0x m) t))
             :qy-expr '(+ q0y (* (e/div p0y m) t))
             :px-expr 'p0x
             :py-expr 'p0y}]
    (is (zero? (core/calculate-analytical-fitness law ds)))))
