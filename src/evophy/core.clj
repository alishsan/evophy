
(ns evophy.core
  (:gen-class)
  (:require [emmy.env :as e]
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


(defn evolve-generation
  [population data population-size]
  (let [elite-n (max 1 (quot population-size 5))
        branch-factor (long (Math/ceil (/ population-size (double elite-n))))]
    (->> population
         (map (fn [ind]
                (assoc ind :fitness (calculate-fitness (:dq ind) (:dp ind) data))))
         (sort-by :fitness >)
         (take elite-n)
         (mapcat (fn [parent]
                   ;; Create a vector containing the parent + mutants
                   (cons parent 
                         (take (dec branch-factor) 
                               (repeatedly #(mutate-individual parent))))))
         (take population-size)
         (vec)))) ;; Ensure it returns a vector for the next generation

(defn -main [& _args]
  (timbre/merge-config! {:min-level :warn})
  (let [data (vec (generate-data 1.0 2.0 1.0 0.0 0.1 100))
        initial (vec (repeatedly 50 random-individual))
        final-pop (reduce (fn [pop _] (evolve-generation pop data 50))
                          initial
                          (range 50))
        best (->> final-pop
                  (map (fn [ind]
                         (assoc ind :fitness (calculate-fitness (:dq ind) (:dp ind) data))))
                  (sort-by :fitness >)
                  first)]
    (println "dq/dt:" (:dq best))
    (println "dp/dt:" (:dp best))
    (println "fitness:" (:fitness best))
    (println "predictions:" (evaluate-predictions (:dq best) (:dp best) data))))
