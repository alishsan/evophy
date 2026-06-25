(ns evophy.mcts
  "MCTS builds analytical genomes; GP selects and mutates each generation.

  Search can stop early (Ctrl+C or [[request-stop!]]) and return the best genome found so far."
  (:require [evophy.core :as core])
  (:import (sun.misc Signal SignalHandler)))

(def ^:private hole '__hole__)
(def ^:private default-exploration 1.41)
(def ^:private max-expr-depth 4)
(def ^:private repair-genetic-samples 6)

(def ^:private analytical-phases [:qx :qy :px :py])
(def ^:private phase->expr-key
  {:qx :qx-expr :qy :qy-expr :px :px-expr :py :py-expr})
(def ^:private expr-key->phase
  (zipmap (vals phase->expr-key) (keys phase->expr-key)))

(defonce ^:private stop-flag (atom false))
(defonce ^:private handler-installed? (atom false))
(defonce ^:private last-search-stats (atom nil))

(defn stop-requested? []
  (or @stop-flag (Thread/interrupted)))

(defn request-stop! []
  (reset! stop-flag true)
  (println "\n[evophy] Stop requested — finishing current step, then returning best-so-far…"))

(defn clear-stop! []
  (reset! stop-flag false)
  (Thread/interrupted) ;; clear interrupt flag if set
  false)

(defn last-search-result []
  @last-search-stats)

(defn install-stop-handler!
  "Register SIGINT (Ctrl+C) to call [[request-stop!]] without killing the JVM immediately."
  []
  (when-not @handler-installed?
    (try
      (Signal/handle "INT"
                     (proxy [SignalHandler] []
                       (handle [_] (request-stop!))))
      (reset! handler-installed? true)
      (catch Throwable _
        nil))))

(defn- hole? [x] (= x hole))

(defn- find-hole-path [expr]
  (cond
    (hole? expr) []
    (sequential? expr)
    (let [args (vec (rest expr))]
      (some (fn [i]
              (when-let [sub (find-hole-path (args i))]
                (cons i sub)))
            (range (count args))))
    :else nil))

(defn- replace-at-path [expr path replacement]
  (if (empty? path)
    replacement
    (let [op (first expr)
          args (vec (rest expr))
          [i & rp] path]
      (cons op (assoc args i (replace-at-path (args i) rp replacement))))))

(defn- depth-at-hole [expr path]
  (loop [d max-expr-depth remaining path tree expr]
    (if (empty? remaining)
      d
      (let [[i & rp] remaining
            args (vec (rest tree))]
        (recur (dec d) rp (args i))))))

(defn- legal-actions [depth vars]
  (let [terminals (for [s (concat vars core/constants)] [:term s])
        ops (if (or (zero? depth) (< (rand) 0.3))
              []
              (concat (for [op core/unary-ops] [:unary op])
                      (for [op core/ops
                            :when (not (contains? core/unary-ops op))]
                        [:binary op])))]
    (if (empty? ops) terminals (vec (concat terminals ops)))))

(defn- legal-actions-deterministic [depth vars]
  (let [terminals (for [s (concat vars core/constants)] [:term s])
        ops (if (zero? depth)
              []
              (concat (for [op core/unary-ops] [:unary op])
                      (for [op core/ops
                            :when (not (contains? core/unary-ops op))]
                        [:binary op])))]
    (if (empty? ops) terminals (vec (concat terminals ops)))))

(defn- apply-action [expr path action]
  (case (first action)
    :term (replace-at-path expr path (second action))
    :unary (replace-at-path expr path (list (second action) hole))
    :binary (replace-at-path expr path (list (second action) hole hole))))

(defn- complete-expr [expr vars]
  (if-let [path (find-hole-path expr)]
    (let [depth (depth-at-hole expr path)
          action (rand-nth (legal-actions depth vars))]
      (complete-expr (apply-action expr path action) vars))
    expr))

(defn- complete-expr-deterministic [expr vars]
  (if-let [path (find-hole-path expr)]
    (let [depth (depth-at-hole expr path)
          actions (legal-actions-deterministic depth vars)]
      (when (seq actions)
        (complete-expr-deterministic (apply-action expr path (first actions)) vars)))
    expr))

(defn- initial-state []
  (into {:phase :qx}
        (map (fn [k] [k hole]) core/analytical-expr-keys)))

(defn- current-phase [state] (:phase state))

(defn- current-expr [state]
  (get state (phase->expr-key (current-phase state))))

(defn- advance-phase [state]
  (if (find-hole-path (current-expr state))
    state
    (let [phase (:phase state)
          idx (.indexOf (clojure.core/vec analytical-phases) phase)]
      (if (and (>= idx 0) (< idx (dec (count analytical-phases))))
        (assoc state :phase (nth analytical-phases (inc idx)))
        (assoc state :phase :done)))))

(defn- terminal? [state]
  (and (= :done (:phase state))
       (every? #(not (find-hole-path (get state %))) core/analytical-expr-keys)))

(defn- update-current [state expr]
  (assoc state (phase->expr-key (current-phase state)) expr))

(defn- state-action-pairs [state]
  (when-let [path (find-hole-path (current-expr state))]
    (let [depth (depth-at-hole (current-expr state) path)
          vars core/analytical-vars]
      (map (fn [action] {:path path :action action})
           (legal-actions depth vars)))))

(defn- apply-state-action [state {:keys [path action]}]
  (-> state
      (update-current (apply-action (current-expr state) path action))
      advance-phase))

(defn- rollout-state [state]
  (loop [s state]
    (if (terminal? s)
      s
      (if-let [pairs (seq (state-action-pairs s))]
        (recur (apply-state-action s (rand-nth pairs)))
        (let [vars core/analytical-vars]
          (recur (into (assoc s :phase :done)
                       (map (fn [k] [k (complete-expr (get s k) vars)])
                            core/analytical-expr-keys))))))))

(defn- state->individual [state]
  (-> (into {:strategy :analytical}
            (map (fn [k] [k (get state k)]) core/analytical-expr-keys))
      (core/ensure-symbol-coverage core/analytical-expr-keys
                                   core/required-analytical-symbols)
      core/ensure-analytical-uses-t))

(defonce ^:private fitness-cache (atom {}))

(defn- fitness-for-individual [ind datasets & {:keys [evaluation aggregate percentile]}]
  (if (core/genome-valid? ind)
    (let [k (core/individual-genome-key ind)]
      (if-let [hit (@fitness-cache k)]
        hit
        (let [fit (core/calculate-fitness-scenarios ind datasets
                                                   :evaluation (or evaluation :de-driven)
                                                   :aggregate (or aggregate :min)
                                                   :percentile (or percentile 0.1))]
          (swap! fitness-cache assoc k fit)
          fit)))
    0.0))

(defn- adversarial-fitness [ind adversarial-dataset all-datasets fit-opts]
  (let [leg (core/first-law-legacy ind)
        adv (when (and (= :analytical (:strategy leg)) adversarial-dataset)
              (core/calculate-analytical-fitness-at-horizon
               leg adversarial-dataset core/repair-horizon-frac))
        full (fitness-for-individual ind all-datasets fit-opts)]
    (if (and adv (pos? adv))
      (+ (* 0.75 adv) (* 0.25 full))
      full)))

(defn- repair-worst-scenario-fitness [ind dataset]
  (let [leg (core/first-law-legacy ind)]
    (if (and dataset (= :analytical (:strategy leg)))
      (core/calculate-analytical-fitness-at-horizon
       leg dataset core/repair-horizon-frac)
      0.0)))

(defn- score-state [state datasets]
  (fitness-for-individual (state->individual (rollout-state state)) datasets))

(defn- ucb1 [child-visits child-value-sum parent-visits c]
  (if (zero? child-visits)
    Double/POSITIVE_INFINITY
    (+ (/ child-value-sum child-visits)
       (* c (Math/sqrt (/ (Math/log (max 1 parent-visits)) child-visits))))))

(defn- select-child [parent-atom c]
  (when-let [children (seq (:children @parent-atom))]
    (->> children
         (map (fn [[action child-atom]]
                [action child-atom
                 (ucb1 (:visits @child-atom) (:value-sum @child-atom)
                       (:visits @parent-atom) c)]))
         (sort-by #(nth % 2) #(compare %2 %1))
         first)))

(defn- backprop! [node-atom score]
  (loop [n node-atom]
    (when n
      (swap! n update :visits inc)
      (swap! n update :value-sum + score)
      (recur (:parent @n)))))

(defn- note-best! [best-atom ind datasets]
  (when (core/genome-valid? ind)
    (let [fit (fitness-for-individual ind datasets)]
      (when (> fit (:fitness @best-atom))
        (swap! best-atom assoc :fitness fit :individual ind)))))

(defn- expand! [node datasets _c best-atom]
  (when-let [pair (first (:untried @node))]
    (swap! node update :untried rest)
    (let [child-state (apply-state-action (:state @node) pair)
          action (:action pair)
          child (atom {:state child-state
                       :visits 0 :value-sum 0.0
                       :children {}
                       :untried (vec (or (state-action-pairs child-state) []))
                       :parent node})]
      (swap! node assoc-in [:children action] child)
      (note-best! best-atom (state->individual (rollout-state child-state)) datasets)
      (let [score (score-state child-state datasets)]
        (backprop! child score)))))

(defn- traverse [root datasets c best-atom]
  (loop [node root]
    (cond
      (seq (:untried @node))
      (do (expand! node datasets c best-atom) node)

      (seq (:children @node))
      (let [[_ child-atom _] (select-child node c)]
        (recur child-atom))

      :else node)))

(defn- resolve-individual [best-atom]
  (or (when-let [i (:individual @best-atom)]
        (when (core/genome-valid? i) i))
      (core/random-analytical-individual)))

(defn search-analytical-individual
  "Run MCTS for up to max-simulations tree traversals (or until [[stop-requested?]]).

  Returns the best analytical genome found. Metadata is stored in [[last-search-result]].
  Optional :should-stop? fn (default [[stop-requested?]]), :on-progress fn called as {:sim n :best-fitness f}."
  [datasets max-simulations & {:keys [should-stop? on-progress]
                               :or {should-stop? stop-requested?}}]
  (when (pos? max-simulations)
    (reset! fitness-cache {})
    (let [c default-exploration
          root-state (initial-state)
          root (atom {:state root-state
                      :visits 0 :value-sum 0.0
                      :children {}
                      :untried (vec (or (state-action-pairs root-state) []))
                      :parent nil})
          best (atom {:fitness -1.0 :individual nil})
          ran (loop [sim 0]
                (if (or (>= sim max-simulations) (should-stop?))
                  sim
                  (do
                    (traverse root datasets c best)
                    (when on-progress
                      (on-progress {:sim (inc sim) :best-fitness (:fitness @best)}))
                    (recur (inc sim)))))]
      (let [stopped? (should-stop?)]
        (when stopped?
          (println (str "[MCTS] stopped after " ran " simulation"
                        (when (not= 1 ran) "s")
                        " (best fitness " (format "%.4f" (:fitness @best)) ")")))
        (let [ind (resolve-individual best)]
          (reset! last-search-stats
                  {:individual ind
                   :simulations-run ran
                   :max-simulations max-simulations
                   :best-fitness (:fitness @best)
                   :stopped? stopped?})
          (if (core/genome-valid? ind) ind (core/random-analytical-individual)))))))

(defn generation-inject
  "Inject up to n MCTS individuals; stops early if [[stop-requested?]]."
  [datasets max-simulations n]
  (loop [i 0 acc []]
    (if (or (>= i n) (stop-requested?))
      acc
      (recur (inc i) (conj acc (search-analytical-individual datasets max-simulations))))))

;; ── Adversarial repair (patch weakest slot on HoF / elite) ───────────────────

(defn- repair-state-from-individual [ind expr-key]
  (let [leg (core/first-law-legacy ind)
        phase (expr-key->phase expr-key)]
    (into {:phase phase :repair-key expr-key}
          (map (fn [k]
                 [k (if (= k expr-key)
                      hole
                      (core/normalize-expr (get leg k)))])
               core/analytical-expr-keys))))

(defn- repair-terminal? [state]
  (not (find-hole-path (get state (:repair-key state)))))

(defn- repair-advance-phase [state]
  (if (find-hole-path (get state (:repair-key state)))
    state
    (assoc state :phase :done)))

(defn- repair-state-action-pairs [state]
  (when-let [path (find-hole-path (get state (:repair-key state)))]
    (let [expr (get state (:repair-key state))
          depth (depth-at-hole expr path)
          vars core/analytical-vars]
      (map (fn [action] {:path path :action action})
           (legal-actions depth vars)))))

(defn- apply-repair-action [state {:keys [path action]}]
  (let [expr-key (:repair-key state)
        expr (get state expr-key)]
    (-> state
        (assoc expr-key (apply-action expr path action))
        repair-advance-phase)))

(defn- repair-finish-state [state]
  (if (repair-terminal? state)
    state
    (let [expr-key (:repair-key state)
          vars core/analytical-vars
          expr (or (complete-expr-deterministic (get state expr-key) vars)
                   (get state expr-key))]
      (assoc state expr-key expr))))

(defn- wrap-repaired-slot-expr [expr]
  (let [ex (core/normalize-expr expr)]
    (if (and (sequential? ex)
             (= 'e/if (first ex))
             (= '(neg? energy) (second ex)))
      ex
      (list 'e/if '(neg? energy) ex ex))))

(defn- repair-state->individual [state base-ind]
  (let [leg (core/first-law-legacy base-ind)
        repair-key (:repair-key state)
        repaired (into leg
                       (map (fn [k]
                              (let [ex (core/normalize-expr (get state k))]
                                [k (if (= k repair-key)
                                     (if core/*both-regimes?*
                                       (wrap-repaired-slot-expr ex)
                                       ex)
                                     ex)]))
                            core/analytical-expr-keys))
        repaired (-> repaired
                     (core/ensure-symbol-coverage core/analytical-expr-keys
                                                  core/required-analytical-symbols)
                     core/ensure-analytical-uses-t)]
    (core/apply-analytical-repair base-ind repaired)))

(defn- repair-genetic-eval [state base-ind adversarial-dataset all-datasets fit-opts]
  "From an expanded MCTS repair node, finish the slot deterministically and score
   the seed plus GP mutants (replaces random hole-filling rollout)."
  (let [finished (repair-finish-state state)
        seed-ind (repair-state->individual finished base-ind)
        seed-leg (core/first-law-legacy seed-ind)
        mutate (fn []
                 (core/apply-analytical-repair base-ind
                                               (core/gp-mutate-analytical seed-leg)))
        candidates (cons seed-ind
                           (repeatedly (dec repair-genetic-samples) mutate))
        scored (keep (fn [ind]
                       (when (core/genome-valid? ind)
                         {:individual ind
                          :score (adversarial-fitness ind adversarial-dataset
                                                      all-datasets fit-opts)}))
                     candidates)]
    (or (apply max-key :score scored)
        {:individual seed-ind
         :score (adversarial-fitness seed-ind adversarial-dataset
                                     all-datasets fit-opts)})))

(defn- note-best-repair! [best-atom ind adversarial-dataset]
  (when (and adversarial-dataset (core/genome-valid? ind))
    (let [fit (repair-worst-scenario-fitness ind adversarial-dataset)]
      (when (> fit (:repair-fitness @best-atom))
        (swap! best-atom assoc :repair-fitness fit :individual ind)))))

(defn- expand-repair! [node adversarial-dataset all-datasets fit-opts best-atom]
  (when-let [pair (first (:untried @node))]
    (swap! node update :untried rest)
    (let [child-state (apply-repair-action (:state @node) pair)
          action (:action pair)
          child (atom {:state child-state
                       :visits 0 :value-sum 0.0
                       :children {}
                       :untried (vec (or (repair-state-action-pairs child-state) []))
                       :parent node})]
      (swap! node assoc-in [:children action] child)
      (let [{:keys [individual score]} (repair-genetic-eval child-state
                                                             (:base-ind (:state @node))
                                                             adversarial-dataset
                                                             all-datasets
                                                             fit-opts)]
        (note-best-repair! best-atom individual adversarial-dataset)
        (backprop! child score)))))

(defn- traverse-repair [root adversarial-dataset all-datasets fit-opts best-atom]
  (loop [node root]
    (cond
      (seq (:untried @node))
      (do (expand-repair! node adversarial-dataset all-datasets fit-opts best-atom) node)

      (seq (:children @node))
      (let [[_ child-atom _] (select-child node default-exploration)]
        (recur child-atom))

      :else node)))

(defn search-repair-individual
  "MCTS adversarial repair: keep HoF/elite fixed except one weak slot, re-search that slot
   against the individual's worst scenario. Returns repaired individual or base-ind on failure.

   Optional :fitness-opts map passed to aggregate fitness; :on-progress {:sim :best-fitness}."
  [base-ind all-datasets max-simulations & {:keys [should-stop? on-progress fitness-opts quiet?]
                                            :or {should-stop? stop-requested?
                                                 fitness-opts {}
                                                 quiet? false}}]
  (if (zero? max-simulations)
    base-ind
    (if-let [target (core/adversarial-repair-target base-ind all-datasets)]
      (let [{:keys [dataset expr-key scenario-id]} target
            fit-opts (merge {:evaluation :de-driven :aggregate :min :percentile 0.1}
                            fitness-opts)
            root-state (assoc (repair-state-from-individual base-ind expr-key)
                              :base-ind base-ind)
            _ (when (and scenario-id (not quiet?))
                (println (format "  [MCTS repair] scenario=%s slot=%s horizon=%.0f%%"
                                 (name scenario-id) (name expr-key)
                                 (* 100.0 core/repair-horizon-frac))))
            parent-repair-fit (repair-worst-scenario-fitness base-ind dataset)
            root (atom {:state root-state
                        :visits 0 :value-sum 0.0
                        :children {}
                        :untried (vec (or (repair-state-action-pairs root-state) []))
                        :parent nil})
            best (atom {:repair-fitness parent-repair-fit :individual base-ind})
            ran (loop [sim 0]
                  (if (or (>= sim max-simulations) (should-stop?))
                    sim
                    (do
                      (traverse-repair root dataset all-datasets fit-opts best)
                      (when on-progress
                        (on-progress {:sim (inc sim)
                                      :best-fitness (:repair-fitness @best)}))
                      (recur (inc sim)))))]
        (let [ind (or (:individual @best) base-ind)
              improved? (> (:repair-fitness @best) parent-repair-fit)
              result (if improved? ind base-ind)
              stopped? (should-stop?)]
          (reset! last-search-stats
                  {:individual result
                   :simulations-run ran
                   :max-simulations max-simulations
                   :best-fitness (:repair-fitness @best)
                   :parent-repair-fitness parent-repair-fit
                   :improved? improved?
                   :stopped? stopped?
                   :repair? true
                   :repair-target target})
          (when (and (not quiet?) stopped? (pos? ran))
            (println (str "[MCTS repair] stopped after " ran " simulation"
                          (when (not= 1 ran) "s")
                          " (repair fitness " (format "%.4f" (:repair-fitness @best))
                          " vs parent " (format "%.4f" parent-repair-fit) ")")))
          (when (and (not quiet?) (pos? ran) (not improved?))
            (println (format "  [MCTS repair] no improvement on %.0f%% horizon (%.4f <= %.4f)"
                             (* 100.0 core/repair-horizon-frac)
                             (:repair-fitness @best) parent-repair-fit)))
          result))
      base-ind)))

(defn generation-repair-inject
  "Up to n adversarial repairs from one source individual (typically HoF).
   Skips repairs that do not improve worst-scenario fitness at [[core/repair-horizon-frac]]."
  [base-ind datasets max-simulations n]
  (loop [i 0 acc []]
    (if (or (>= i n) (stop-requested?))
      acc
      (let [repaired (search-repair-individual base-ind datasets max-simulations)
            stats (last-search-result)]
        (recur (inc i)
               (if (:improved? stats)
                 (conj acc repaired)
                 acc))))))
