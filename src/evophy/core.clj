(ns evophy.core
  (:require [emmy.env :as e]
            [taoensso.timbre :as timbre]))

(defn ho-hamiltonian [m k]
  (fn [[_ q p]]
    ;; Use emmy.generic arithmetic so H works under automatic differentiation
    ;; inside [[Hamiltonian->state-derivative]] (Differential tags on q, p).
    (e/+ (e// (e/square p) (e/* 2 m))
         (e/* (e// 1 2) k (e/square q)))))

(defn ho-state-derivative [m k]
  (e/Hamiltonian->state-derivative (ho-hamiltonian m k)))

(defn generate-data [m k q0 p0 dt steps]
  (let [h          (ho-hamiltonian m k)
        initial    (e/->H-state 0 q0 p0)
        t-final    (* dt steps)
        ;; derivative-args are coerced with double-array in Emmy ODE—use [m k], not [h].
        trajectory (e/integrate-state-derivative ho-state-derivative [m k] initial t-final dt)]
    (map (fn [state]
           (let [[t q p] (vec state)]
             {:t t :q q :p p :energy (h state)}))
         trajectory)))

(defn candidate-law [q p]
  ;; This represents a 'guess' at the physics law.
  ;; A true AlphaZero would generate these via search.
  ;; For now, let's test a simple guess: p^2 + q^2
  (+ (* 0.5 (e/square p))
     (* 0.5 (e/square q))))

(defn evaluate-law [law-fn data]
  (let [results (map (fn [{:keys [q p]}] (law-fn q p)) data)
        mean    (/ (reduce + results) (count results))
        ;; Variance measures how 'constant' the law is.
        ;; In physics, zero variance = a discovered law!
        variance (/ (reduce + (map #(e/square (- % mean)) results))
                    (count results))]
    {:mean mean :variance variance}))

(def ops '[+ - * e/square])
(def vars '[q p])
(def constants '[0.5 1.0 2.0])

(defn random-atom []
  (rand-nth (concat vars constants)))

(defn random-expression [depth]
  (if (or (zero? depth) (< (rand) 0.3))
    (random-atom)
    (let [op (rand-nth ops)]
      (if (= op 'e/square)
        (list op (random-expression (dec depth)))
        (list op (random-expression (dec depth)) (random-expression (dec depth)))))))

(defn calculate-fitness [expr data]
  (try
    (let [;; Compile the S-expression into a Clojure function
          func (eval `(fn [~'q ~'p] ~expr))
          values (map (fn [{:keys [q p]}] (func q p)) data)
          mean (/ (reduce + values) (count values))
          variance (/ (reduce + (map #(e/square (- % mean)) values)) (count values))
          cv (/ (Math/sqrt variance) (+ (Math/abs mean) 1e-6))]
      ;; Higher fitness when coefficient of variation is low (stable law vs scale).
      (/ 1.0 (+ cv 1e-6)))
    (catch Exception e 0))) ;; If the formula does something illegal (div by zero), fitness is 0

(defn mutate [expr]
  (if (< (rand) 0.2)
    (random-expression 2) ;; 20% chance to totally replace a branch
    (if (coll? expr)
      (map #(if (coll? %) (mutate %) %) expr)
      expr)))

(defn evolve-generation
  ([population data]
   (evolve-generation population data 30))
  ([population data population-size]
   (let [elite-n (max 1 (quot population-size 5))
         branch-factor (long (Math/ceil (/ population-size (double elite-n))))]
     (->> population
          (map (fn [expr] {:expr expr :fitness (calculate-fitness expr data)}))
          (sort-by :fitness >)
          (take elite-n)
          (mapcat (fn [parent]
                    (cons (:expr parent)
                          (take (dec branch-factor)
                                (repeatedly #(mutate (:expr parent)))))))
          (take population-size)))))

(defn -main [& args]
  ;; Emmy logs ODE compile at INFO; keep default for REPL, quiet for `lein run`.
  (timbre/merge-config! {:min-level :warn})
  (let [data      (generate-data 1.0 1.0 1.0 0.0 0.1 100)
        initial   (vec (repeatedly 50 #(random-expression 4)))
        final-pop (reduce (fn [pop _] (evolve-generation pop data 50))
                          initial
                          (range 50))
        best      (->> final-pop
                       (map (fn [expr] {:expr expr :fitness (calculate-fitness expr data)}))
                       (sort-by :fitness >)
                       first)]
    (println (:expr best))))
