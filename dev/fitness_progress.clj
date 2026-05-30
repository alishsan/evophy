(ns fitness-progress
  (:require [evophy.core :as c]))

;; ---------------------------------------------------------------------------
;; Helpers
;; ---------------------------------------------------------------------------

(defn- sparkline
  "Map a sequence of doubles to an 8-bar Unicode sparkline string."
  [values]
  (let [bars "▁▂▃▄▅▆▇█"
        lo   (apply min values)
        hi   (apply max values)
        span (- hi lo)]
    (apply str
           (map (fn [v]
                  (let [idx (if (zero? span)
                              4
                              (int (* 7.0 (/ (- v lo) span))))]
                    (nth bars (max 0 (min 7 idx)))))
                values))))

(defn- pad-right [s n]
  (let [s (str s)]
    (if (>= (count s) n) s (str s (apply str (repeat (- n (count s)) " "))))))

(defn- fmt  [x] (if x (format "%.4f" (double x)) "  n/a "))
(defn- fmt+ [x] (if (and x (pos? x)) (format "%+.4f" (double x)) "   n/a"))

;; ---------------------------------------------------------------------------
;; Table renderer
;; ---------------------------------------------------------------------------

(defn- print-table [history]
  (let [has-hof?  (some :hof-best history)
        has-eval? (some :eval-best history)
        header    (str (pad-right "gen" 6)
                       (pad-right "total-gen" 10)
                       (when has-hof?  (pad-right "hof-best†" 11))
                       (when has-eval? (pad-right "eval-best*" 12))
                       (pad-right "train-best" 12)
                       (pad-right "mean" 10))]
    (println header)
    (println (apply str (repeat (count header) "-")))
    (doseq [{:keys [gen total-gen best mean eval-best hof-best]} history]
      (println (str (pad-right gen 6)
                    (pad-right total-gen 10)
                    (when has-hof?  (pad-right (fmt hof-best) 11))
                    (when has-eval? (pad-right (fmt eval-best) 12))
                    (pad-right (fmt best) 12)
                    (pad-right (fmt mean) 10))))
    (println)
    (when has-hof?
      (println "  † hof-best  = hall-of-fame: best individual ever seen (should only increase)"))
    (when has-eval?
      (println "  * eval-best = best elite individual this generation on fixed scenarios (noisy)"))))

;; ---------------------------------------------------------------------------
;; Sparkline charts
;; ---------------------------------------------------------------------------

(defn- print-chart [history]
  (let [has-hof?  (some :hof-best history)
        has-eval? (some :eval-best history)
        n         (count history)]
    (when has-hof?
      (let [hofs (mapv #(or (:hof-best %) 0.0) history)
            hi   (apply max hofs)
            lo   (apply min hofs)
            last-entry (last history)
            first-entry (first history)
            delta (when (and (:hof-best first-entry) (:hof-best last-entry))
                    (- (:hof-best last-entry) (:hof-best first-entry)))]
        (println (str "\nHall-of-fame best (fixed scenarios) — should be monotonically non-decreasing"))
        (println (str "  hi=" (fmt hi) "  lo=" (fmt lo)
                      (when delta (str "  total-gain=" (fmt+ delta)))))
        (println (str "  " (sparkline hofs)))))
    (when has-eval?
      (let [evals (mapv #(or (:eval-best %) 0.0) history)
            hi    (apply max evals)
            lo    (apply min evals)]
        (println (str "\nEval-best per generation (fixed scenarios, noisy)"))
        (println (str "  hi=" (fmt hi) "  lo=" (fmt lo)))
        (println (str "  " (sparkline evals)))))
    (let [bests (mapv :best history)
          means (mapv :mean history)
          hi    (apply max bests)
          lo    (apply min bests)]
      (println (str "\nTrain-best (random scenarios, not comparable across gens)"))
      (println (str "  hi=" (fmt hi) "  lo=" (fmt lo)))
      (println (str "  " (sparkline bests)))
      (println "\nMean train fitness")
      (println (str "  " (sparkline means))))))

;; ---------------------------------------------------------------------------
;; Summary
;; ---------------------------------------------------------------------------

(defn- print-summary [history]
  (let [last-entry  (last history)
        first-entry (first history)
        has-hof?    (some :hof-best history)
        hof-delta   (when (and has-hof? (:hof-best first-entry) (:hof-best last-entry))
                      (- (:hof-best last-entry) (:hof-best first-entry)))]
    (println "\nSummary")
    (println (str "  generations logged : " (count history)))
    (when has-hof?
      (println (str "  hall-of-fame best  : " (fmt (:hof-best last-entry)))))
    (when hof-delta
      (println (str "  improvement (hof)  : " (fmt+ hof-delta))))
    (println (str "  eval-best (last)   : " (fmt (:eval-best last-entry))))
    (println (str "  train-best (last)  : " (fmt (:best last-entry))))
    (println (str "  mean (last)        : " (fmt (:mean last-entry))))))

;; ---------------------------------------------------------------------------
;; Entry point
;; ---------------------------------------------------------------------------

(let [pop-path  "data/population.edn"
      hist-path (c/history-path-for pop-path)
      loaded    (c/load-history hist-path)]
  (if (or (nil? loaded) (empty? loaded))
    (println (str "No fitness history at " hist-path " — run lein run first."))
    (do
      (print-table loaded)
      (print-chart loaded)
      (print-summary loaded))))
