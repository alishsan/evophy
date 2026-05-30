(require '[evophy.core :as c])

(defn- pct [x] (format "%.0f%%" (* 100.0 (double x))))

(defn- print-motif [{:keys [slot signature share avg-fitness tags gloss motif]}]
  (println (str "  [" (name slot) "] share=" (pct share)
                " avg-fit=" (format "%.3f" (double avg-fitness))
                (when (seq tags) (str " tags=" tags))
                (when gloss (str "  → " gloss))))
  (println (str "       " signature)))

(defn- print-report [report]
  (println "Population motif report")
  (println "  individuals:" (:n-population report)
             "  elite mined:" (:n-elite report)
             "  elite min fitness:" (format "%.4f" (double (:elite-min-fitness report))))
  (println)
  (if (empty? (:physics-hits report))
    (println "Physics-shaped motifs (exact Hamilton matches): none")
    (do
      (println "Physics-shaped motifs (exact Hamilton matches):")
      (doseq [m (:physics-hits report)] (print-motif m))))
  (println)
  (println "Top recurring micro-heuristics (by elite share × fitness):")
  (if (empty? (:motifs report))
    (println "  (none above min-share threshold)")
    (doseq [m (:motifs report)] (print-motif m))))

(let [path "data/population.edn"
      loaded (c/load-population path)]
  (if-not loaded
    (println "No checkpoint at" path "- run lein run first")
    (print-report (c/population->motif-report (:population loaded)
                                              :elite-frac 0.5
                                              :min-share 0.15
                                              :top 20))))
