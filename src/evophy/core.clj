
(ns evophy.core
  (:gen-class)
  (:require [clojure.edn :as edn]
            [emmy.env :as e]
            [taoensso.timbre :as timbre]))

(defn ho-hamiltonian [m k]
  (fn [[_ q p]]
    (e/+ (e// (e/square p) (e/* 2 m))
         (e/* (e// 1 2) k (e/square q)))))

(defn ho-state-derivative [m k]
  (e/Hamiltonian->state-derivative (ho-hamiltonian m k)))

(defn generate-data [m k q0 p0 dt steps]
  (let [h (ho-hamiltonian m k)
        initial (e/->H-state 0 q0 p0)
        t-final (* dt steps)
        trajectory (e/integrate-state-derivative ho-state-derivative [m k] initial t-final dt)]
    (map (fn [state]
           (let [[t q p] (vec state)]
             {:t t :q q :p p :energy (h state)}))
         trajectory)))

(defn scenario-data
  "Trajectory for one training scenario (map with :m :k :q0 :p0 :dt :steps)."
  [{:keys [m k q0 p0 dt steps]}]
  (vec (generate-data m k q0 p0 dt steps)))

(def default-scenarios
  "Diverse harmonic-oscillator setups—like many AlphaZero games with different openings."
  [{:id :displaced-q :m 1.0 :k 1.0 :q0 1.0 :p0 0.0 :dt 0.1 :steps 100}
   {:id :displaced-p :m 1.0 :k 1.0 :q0 0.0 :p0 1.0 :dt 0.1 :steps 100}
   {:id :heavy-m      :m 2.0 :k 1.0 :q0 0.8 :p0 0.5 :dt 0.1 :steps 100}
   {:id :stiff-k      :m 1.0 :k 2.0 :q0 0.5 :p0 0.5 :dt 0.1 :steps 100}
   {:id :mixed-phase  :m 1.0 :k 1.0 :q0 0.7 :p0 0.7 :dt 0.05 :steps 80}])

(defn scenarios->datasets
  [scenarios]
  (mapv scenario-data scenarios))

(defn- compile-rate-fn [expr]
  (binding [*ns* (the-ns 'evophy.core)]
    (eval `(fn [~'q ~'p] ~expr))))

(defn- predict-step [dq-fn dp-fn {:keys [q p]} dt]
  (let [dq (double (dq-fn q p))
        dp (double (dp-fn q p))]
    {:q (+ (double q) (* dq dt))
     :p (+ (double p) (* dp dt))}))

(defn- prediction-errors [dq-fn dp-fn data]
  (reduce
   (fn [{:keys [sq-q sq-p n] :as acc} [s0 s1]]
     (let [dt (- (:t s1) (:t s0))
           pred (predict-step dq-fn dp-fn s0 dt)]
       (-> acc
           (update :sq-q + (double (e/square (- (:q pred) (:q s1)))))
           (update :sq-p + (double (e/square (- (:p pred) (:p s1)))))
           (update :n inc))))
   {:sq-q 0.0 :sq-p 0.0 :n 0}
   (partition 2 1 data)))

(defn evaluate-predictions
  "One-step Euler forecast error along consecutive trajectory samples."
  [dq-expr dp-expr data]
  (let [dq-fn (compile-rate-fn dq-expr)
        dp-fn (compile-rate-fn dp-expr)
        {:keys [sq-q sq-p n]} (prediction-errors dq-fn dp-fn data)]
    (if (zero? n)
      {:n-steps 0 :mse-q Double/POSITIVE_INFINITY :mse-p Double/POSITIVE_INFINITY
       :mse Double/POSITIVE_INFINITY}
      (let [n (double n)]
        {:n-steps (long n)
         :mse-q (/ sq-q n)
         :mse-p (/ sq-p n)
         :mse (/ (+ sq-q sq-p) n)}))))

(def ops '[+ - * e/square])
(def vars '[q p])
(def constants '[0.5 1.0 2.0])

(def operators
  {'+      {:arity 2 :fn e/+}
   '-      {:arity 2 :fn e/-}
   '*      {:arity 2 :fn e/*}
   'e/div  {:arity 2 :fn (fn [a b] (if (< (Math/abs (double b)) 1e-6) 1.0 (e// a b)))}
   'e/sin  {:arity 1 :fn e/sin}
   'e/cos  {:arity 1 :fn e/cos}
   'e/exp  {:arity 1 :fn (fn [x] (let [val (double x)]
                                   (if (> val 20.0) (e/exp 20.0) (e/exp val))))} ;; Cap to avoid Inf
   'e/square {:arity 1 :fn e/square}})

(defn random-atom []
  (rand-nth (concat vars constants)))

(defn random-expression [depth]
  (if (or (zero? depth) (< (rand) 0.3))
    (random-atom)
    (let [op (rand-nth ops)]
      (if (= op 'e/square)
        (list op (random-expression (dec depth)))
        (list op
              (random-expression (dec depth))
              (random-expression (dec depth)))))))

(defn random-individual []
  {:dq (random-expression 4)
   :dp (random-expression 4)})


(def terminals ['q 'p (fn [] (rand-nth [0.5 1.0 2.0]))])

(defn random-tree [max-depth]
  (if (or (zero? max-depth) (< (rand) 0.3))
    ;; Terminal node (Leaf)
    (let [t (rand-nth terminals)]
      (if (fn? t) (t) t))
    ;; Operator node (Branch)
    (let [op-name (rand-nth (keys operators))
          arity   (get-in operators [op-name :arity])
          children (repeatedly arity #(random-tree (dec max-depth)))]
      (cons op-name children))))

(defn calculate-fitness
  "Higher when dq/dt and dp/dt laws forecast the next (q,p) accurately and quickly."
  [dq-expr dp-expr data]
  (try
    (let [dq-fn (compile-rate-fn dq-expr)
          dp-fn (compile-rate-fn dp-expr)
          {:keys [sq-q sq-p n]} (prediction-errors dq-fn dp-fn data)]
      (if (zero? n)
        0
        (let [n (double n)
              rmse (Math/sqrt (/ (+ sq-q sq-p) n))
              quality (/ 1.0 (+ rmse 1e-6))
              pairs (partition 2 1 data)
              t0 (System/nanoTime)]
          (doseq [[s0 s1] pairs]
            (predict-step dq-fn dp-fn s0 (- (:t s1) (:t s0))))
          (let [elapsed-ms (/ (double (- (System/nanoTime) t0)) 1e6)
                speed-factor (/ 1.0 (+ elapsed-ms 0.01))]
            (* quality speed-factor)))))
    (catch Exception _ 0)))

(defn calculate-fitness-scenarios
  "Mean fitness across many trajectories—laws must generalize, not overfit one IC."
  [dq-expr dp-expr datasets]
  (let [fits (mapv #(calculate-fitness dq-expr dp-expr %) datasets)
        n (count fits)]
    (if (zero? n)
      0
      (/ (reduce + 0.0 fits) n))))

(defn mutate [expr]
  (if (< (rand) 0.2)
    (random-expression 2)
    (if (coll? expr)
      (map #(if (coll? %) (mutate %) %) expr)
      expr)))

(defn mutate-individual [{:keys [dq dp]}]
  (if (< (rand) 0.2)
    (random-individual)
    {:dq (mutate dq) :dp (mutate dp)}))


(def default-population-file "data/population.edn")
(def checkpoint-version 1)

(defn- individual-for-save [{:keys [dq dp]}]
  {:dq dq :dp dp})

(defn save-population!
  "Write population to path (creates parent dirs). Omits :fitness; it is recomputed on load."
  [path population & {:keys [generations-run population-size]}]
  (let [f (java.io.File. path)]
    (.mkdirs (.getParentFile f))
    (spit path
          (pr-str {:version checkpoint-version
                   :population-size (or population-size (count population))
                   :generations-run (long (or generations-run 0))
                   :population (mapv individual-for-save population)}))))

(defn load-population
  "Load checkpoint from path, or nil if missing / invalid."
  [path]
  (let [f (java.io.File. path)]
    (when (.exists f)
      (try
        (let [{:keys [version population generations-run population-size]}
              (edn/read-string (slurp path))]
          (when (= version checkpoint-version)
            {:population (vec population)
             :generations-run (long (or generations-run 0))
             :population-size population-size}))
        (catch Exception e
          (timbre/warn "Could not load population from" path "-" (.getMessage e))
          nil)))))

(defn normalize-population-size
  "Pad with random individuals or trim so size matches the configured target."
  [population target-size]
  (let [n (count population)]
    (cond
      (zero? target-size) []
      (<= n target-size) (into population (repeatedly (- target-size n) random-individual))
      :else (subvec (vec population) 0 target-size))))

(defn resolve-initial-population
  [{:keys [fresh? path population-size]}]
  (if fresh?
    {:population (vec (repeatedly population-size random-individual))
     :generations-run 0
     :resumed? false}
    (if-let [{:keys [population generations-run]} (load-population path)]
      {:population population
       :generations-run generations-run
       :resumed? true}
      {:population (vec (repeatedly population-size random-individual))
       :generations-run 0
       :resumed? false})))

(defn parse-args
  [args]
  (loop [opts {:fresh? false
               :path default-population-file
               :generations 50
               :population-size 50}
         xs args]
    (if (empty? xs)
      opts
      (let [[a & more] xs]
        (case a
          "--fresh" (recur (assoc opts :fresh? true) more)
          "--population" (recur (assoc opts :path (first more)) (rest more))
          "--generations" (recur (assoc opts :generations (Long/parseLong (first more))) (rest more))
          "--population-size" (recur (assoc opts :population-size (Long/parseLong (first more))) (rest more))
          (throw (ex-info "Unknown argument (try --fresh --population PATH --generations N --population-size N)"
                          {:arg a})))))))

(defn evolve-generation
  [population datasets population-size]
  (let [elite-n (max 1 (quot population-size 5))
        branch-factor (long (Math/ceil (/ population-size (double elite-n))))]
    (->> population
         (map (fn [ind]
                (assoc ind :fitness (calculate-fitness-scenarios (:dq ind) (:dp ind) datasets))))
         (sort-by :fitness >)
         (take elite-n)
         (mapcat (fn [parent]
                   (cons parent
                         (take (dec branch-factor)
                               (repeatedly #(mutate-individual parent))))))
         (take population-size)
         (vec))))

(defn -main [& args]
  (timbre/merge-config! {:min-level :warn})
  (let [{:keys [fresh? path generations population-size]} (parse-args args)
        scenarios default-scenarios
        datasets (scenarios->datasets scenarios)
        {:keys [population generations-run resumed?]}
        (resolve-initial-population {:fresh? fresh?
                                     :path path
                                     :population-size population-size})
        initial (normalize-population-size population population-size)
        final-pop (reduce (fn [pop _] (evolve-generation pop datasets population-size))
                          initial
                          (range generations))
        total-generations (+ generations-run generations)
        best (->> final-pop
                  (map (fn [ind]
                         (assoc ind :fitness (calculate-fitness-scenarios (:dq ind) (:dp ind) datasets))))
                  (sort-by :fitness >)
                  first)]
    (save-population! path final-pop
                      :generations-run total-generations
                      :population-size population-size)
    (println (if resumed? "resumed from" "started new; saved to") path)
    (println "total generations (cumulative):" total-generations)
    (println "scenarios:" (count scenarios))
    (println "dq/dt:" (:dq best))
    (println "dp/dt:" (:dp best))
    (println "mean fitness:" (:fitness best))
    (doseq [scenario scenarios]
      (let [data (scenario-data scenario)
            metrics (evaluate-predictions (:dq best) (:dp best) data)]
        (println (str "  " (:id scenario) " mse=" (format "%.6f" (:mse metrics))))))))
