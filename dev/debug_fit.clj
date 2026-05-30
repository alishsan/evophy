(ns debug-fit
  (:require [evophy.core :as c]))

(def r-cubed
  '(* (e/sqrt (+ (e/square qx) (e/square qy)))
      (e/square (e/sqrt (+ (e/square qx) (e/square qy))))))

(def true-grav-rates
  {:strategy :differential
   :dqx-expr '(e/div px m)
   :dqy-expr '(e/div py m)
   :dpx-expr (list '* -1.0 (list '* 'alpha (list 'e/div 'qx r-cubed)))
   :dpy-expr (list '* -1.0 (list '* 'alpha (list 'e/div 'qy r-cubed)))})

(println "True rates valid?:" (c/genome-valid? true-grav-rates))

(let [datasets (c/scenarios->datasets c/default-scenarios)]
  (doseq [ds datasets]
    (println "Scenario:" (:id ds) "Fitness:" (c/calculate-fitness true-grav-rates ds)))
  (println "Overall worst-scenario fitness:" (c/calculate-fitness-scenarios true-grav-rates datasets)))
