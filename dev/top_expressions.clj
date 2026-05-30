(ns top-expressions
  (:require [evophy.core :as c]))

(def pop-path      "data/population.edn")
(def top-n         5)
(def eval-timeout-ms 1500)

(defn- safe-fitness
  "Score one individual with a timeout; returns 0.0 on timeout or error."
  [ind datasets]
  (let [f (future
            (try (c/calculate-fitness-scenarios ind datasets)
                 (catch Exception _ 0.0)))]
    (deref f eval-timeout-ms 0.0)))

(defn- print-individual [i ind datasets]
  (try
    (let [eq (c/individual->equations ind :format :plain)]
      (println (str "#" (inc i)
                    "  strategy=" (name (:strategy eq))
                    "  mse=" (format "%.6f" (:mse (:metrics eq)))
                    "  fitness=" (format "%.4f" (:fitness ind))))
      (doseq [[_ line] (sort-by key (:equations eq))]
        (println (str "    " line)))
      (println))
    (catch Exception ex
      (println (str "#" (inc i) "  [display error: " (.getMessage ex) "]"))
      (println))))

(defn- print-section [label inds datasets]
  (println (str "── " label " ──"))
  (if (empty? inds)
    (println "  (none)\n")
    (doseq [[i ind] (map-indexed vector inds)]
      (print-individual i ind datasets))))

(let [loaded (c/load-population pop-path)]
  (if-not loaded
    (println (str "No checkpoint at " pop-path " — run lein run first."))
    (let [datasets (c/scenarios->datasets c/default-scenarios)
          pop      (:population loaded)
          total    (count pop)
          _        (println (str "Scoring " total " individuals in parallel (timeout "
                                 eval-timeout-ms "ms each)..."))
          ;; Parallel scoring; skip genome-invalid individuals immediately.
          ranked   (->> pop
                        (pmap (fn [ind]
                                (assoc ind :fitness
                                       (if (c/genome-valid? ind)
                                         (safe-fitness ind datasets)
                                         0.0))))
                        (sort-by :fitness #(compare %2 %1))
                        vec)
          _        (println (str "Done. Best fitness: "
                                 (format "%.4f" (:fitness (first ranked))) "\n"))
          ;; Deduplicate by genome (same expressions → same fitness → equal maps).
          ;; Then take top-N per strategy from the de-duped list.
          deduped  (distinct ranked)
          top-all  (take top-n deduped)
          top-diff (take top-n (filter #(= :differential (:strategy %)) deduped))
          top-anal (take top-n (filter #(= :analytical   (:strategy %)) deduped))]
      (println (str "Loaded " total
                    " individuals  |  total generations: " (:generations-run loaded) "\n"))
      (print-section (str "Top " top-n " overall (any strategy)") top-all  datasets)
      (print-section (str "Top " top-n " differential")           top-diff datasets)
      (print-section (str "Top " top-n " analytical")             top-anal datasets))))
