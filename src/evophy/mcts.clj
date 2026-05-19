(ns evophy.mcts
  "MCTS builds analytical genomes; GP selects and mutates each generation."
  (:require [evophy.core :as core]))

(def ^:private hole '__hole__)
(def ^:private default-exploration 1.41)
(def ^:private max-expr-depth 4)

(def ^:private analytical-phases [:qx :qy :px :py])
(def ^:private phase->expr-key
  {:qx :qx-expr :qy :qy-expr :px :px-expr :py :py-expr})

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
  (into {:strategy :analytical}
        (map (fn [k] [k (get state k)]) core/analytical-expr-keys)))

(defonce ^:private fitness-cache (atom {}))

(defn- fitness-for-individual [ind datasets]
  (if (core/genome-valid? ind)
    (let [k (core/individual-genome-key ind)]
      (if-let [hit (@fitness-cache k)]
        hit
        (let [fit (core/calculate-fitness-scenarios ind datasets)]
          (swap! fitness-cache assoc k fit)
          fit)))
    0.0))

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
         (sort-by #(nth % 2) >)
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
        (reset! best-atom {:fitness fit :individual ind})))))

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

(defn search-analytical-individual
  [datasets simulations]
  (when (pos? simulations)
    (reset! fitness-cache {})
    (let [c default-exploration
          root-state (initial-state)
          root (atom {:state root-state
                      :visits 0 :value-sum 0.0
                      :children {}
                      :untried (vec (or (state-action-pairs root-state) []))
                      :parent nil})
          best (atom {:fitness -1.0 :individual nil})]
      (dotimes [_ simulations]
        (traverse root datasets c best))
      (or (:individual @best) (core/random-analytical-individual)))))

(defn generation-inject
  [datasets simulations n]
  (vec (repeatedly n #(search-analytical-individual datasets simulations))))
