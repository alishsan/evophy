(ns debug-fit
  (:require [evophy.core :as c]))

;; ── exact seed (r3 derived var) ──────────────────────────────────────────────
(def exact-seed
  {:strategy :differential
   :dqx-expr '(e/div px m)
   :dqy-expr '(e/div py m)
   :dpx-expr '(e/div (* (* -1.0 alpha) qx) r3)
   :dpy-expr '(e/div (* (* -1.0 alpha) qy) r3)})

;; ── exact seed WITHOUT r3 (build r3 explicitly using sqrt)  ──────────────────
(def r-cubed
  '(* (e/sqrt (+ (e/square qx) (e/square qy)))
      (e/square (e/sqrt (+ (e/square qx) (e/square qy))))))

(def true-grav-rates
  {:strategy :differential
   :dqx-expr '(e/div px m)
   :dqy-expr '(e/div py m)
   :dpx-expr (list '* -1.0 (list '* 'alpha (list 'e/div 'qx r-cubed)))
   :dpy-expr (list '* -1.0 (list '* 'alpha (list 'e/div 'qy r-cubed)))})

(def energy-seed
  {:strategy :conserved
   :c-expr '(- (e/div (+ (* px px) (* py py)) (* 2.0 m)) (e/div alpha r))})

(def angular-momentum-seed
  {:strategy :conserved
   :c-expr '(- (* qx py) (* qy px))})

(def analytical-euler-seed
  {:strategy :analytical
   :qx-expr '(+ q0x (* (e/div p0x m) t))
   :qy-expr '(+ q0y (* (e/div p0y m) t))
   :px-expr '(+ p0x (* (e/div (* (* -1.0 alpha) q0x)
                              (* (e/sqrt (+ (* q0x q0x) (* q0y q0y)))
                                 (+ (* q0x q0x) (* q0y q0y))))
                       t))
   :py-expr '(+ p0y (* (e/div (* (* -1.0 alpha) q0y)
                              (* (e/sqrt (+ (* q0x q0x) (* q0y q0y)))
                                 (+ (* q0x q0x) (* q0y q0y))))
                       t))})

(println "=== Seed validity checks ===")
(println "analytical-euler-seed valid?:  " (c/genome-valid? analytical-euler-seed))
(println "exact-seed (r3) valid?:        " (c/genome-valid? exact-seed))
(println "true-grav-rates valid?:        " (c/genome-valid? true-grav-rates))
(println "energy-seed valid?:            " (c/genome-valid? energy-seed))
(println "angular-momentum-seed valid?:  " (c/genome-valid? angular-momentum-seed))
(doseq [[i s] (map-indexed vector c/physics-seeds)]
  (println (str "physics-seeds[" i "] valid?: " (c/genome-valid? s))))

(println)
(println "=== Per-scenario fitness ===")
(let [datasets (c/scenarios->datasets c/default-scenarios)]
  (println "--- exact-seed (r3 derived var) ---")
  (doseq [ds datasets]
    (let [f (c/calculate-fitness exact-seed ds)]
      (println (format "  %-20s fitness=%.4f" (:id ds) f))))
  (println (format "  %-20s fitness=%.4f (worst-scenario agg)"
                   "OVERALL" (c/calculate-fitness-scenarios exact-seed datasets)))

  (println)
  (println "--- analytical first-order Taylor (short horizon 15%) ---")
  (doseq [ds datasets]
    (let [f (c/calculate-fitness analytical-euler-seed ds)]
      (println (format "  %-20s fitness=%.4f" (:id ds) f))))
  (println (format "  %-20s fitness=%.4f (worst-scenario agg)"
                   "OVERALL" (c/calculate-fitness-scenarios analytical-euler-seed datasets)))

  (println)
  (println "--- energy H = (px²+py²)/(2m) - alpha/r ---")
  (doseq [ds datasets]
    (let [f (c/calculate-fitness energy-seed ds)]
      (println (format "  %-20s fitness=%.4f" (:id ds) f))))
  (println (format "  %-20s fitness=%.4f (worst-scenario agg)"
                   "OVERALL" (c/calculate-fitness-scenarios energy-seed datasets)))

  (println)
  (println "--- angular momentum L = qx·py - qy·px ---")
  (doseq [ds datasets]
    (let [f (c/calculate-fitness angular-momentum-seed ds)]
      (println (format "  %-20s fitness=%.4f" (:id ds) f))))
  (println (format "  %-20s fitness=%.4f (worst-scenario agg)"
                   "OVERALL" (c/calculate-fitness-scenarios angular-momentum-seed datasets)))

  (println)
  (println "--- all physics-seeds ---")
  (doseq [[i s] (map-indexed vector c/physics-seeds)]
    (let [f (c/calculate-fitness-scenarios s datasets)]
      (println (format "  seed[%d] strategy=%-12s  fitness=%.4f  valid=%s"
                       i (name (:strategy s)) f (c/genome-valid? s))))))
