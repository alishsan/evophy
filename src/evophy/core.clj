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

(defn- state-sensitivity
  "Mean max(|∂f/∂q|, |∂f/∂p|) approximated with tiny (q,p) bumps; ~0 for literals in (q,p)."
  [func data]
  (let [samples (take 32 data)
        n (count samples)]
    (if (zero? n)
      0.0
      (/ (double
          (reduce
           (fn [acc {:keys [q p]}]
             (+ acc
                (let [q (double q)
                      p (double p)
                      eps (* 1e-7 (max 1.0 (+ (Math/abs q) (Math/abs p))))
                      f0 (double (func q p))
                      dq (/ (Math/abs (- (double (func (+ q eps) p)) f0)) eps)
                      dp (/ (Math/abs (- (double (func q (+ p eps))) f0)) eps)]
                  (max dq dp))))
           0.0
           samples))
         n))))

(defn calculate-fitness [expr data]
  (try
    (let [func (binding [*ns* (the-ns 'evophy.core)]
                  (eval `(fn [~'q ~'p] ~expr)))
          values (map (fn [{:keys [q p]}] (func q p)) data)
          n (count values)]
      (if (zero? n)
        0
        (let [mean (/ (reduce + values) n)
              variance (/ (reduce + (map #(e/square (- % mean)) values)) n)
              cv (/ (Math/sqrt variance) (+ (Math/abs mean) 1e-6))
              sens (state-sensitivity func data)]
          ;; Low CV favors conserved-like quantities; sensitivity zeros out pure
          ;; constants (flat in q,p) without requiring particular symbols in the tree.
          (if (< sens 1e-9)
            0
            (* (/ 1.0 (+ cv 1e-6)) (Math/sqrt sens))))))
    (catch Exception _ 0)))

(defn mutate [expr]
  (if (< (rand) 0.2)
    (random-expression 2) ;; 20% chance to totally replace a branch
    (if (coll? expr)
      (map #(if (coll? %) (mutate %) %) expr)
      expr)))

(defn evolve-generation
  ([population data]
   (evolve-generation population data 50))
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

(defn- take-distinct-by
  "Walk coll in order; keep at most n items whose (key-fn x) has not been seen."
  [key-fn n coll]
  (loop [seen #{} out [] xs (seq coll)]
    (if (or (= (count out) n) (nil? xs))
      out
      (let [x (first xs)
            k (key-fn x)]
        (if (contains? seen k)
          (recur seen out (rest xs))
          (recur (conj seen k) (conj out x) (rest xs)))))))

(defn -main [& _args]
  ;; Emmy logs ODE compile at INFO; keep default for REPL, quiet for `lein run`.
  (timbre/merge-config! {:min-level :warn})
  (let [data      (generate-data 1.0 1.0 1.0 0.0 0.1 100)
        initial   (vec (repeatedly 50 #(random-expression 4)))
        final-pop (reduce (fn [pop _] (evolve-generation pop data 50))
                          initial
                          (range 50))
        top3      (->> final-pop
                       (map (fn [expr] {:expr expr :fitness (calculate-fitness expr data)}))
                       (sort-by :fitness >)
                       (take-distinct-by :expr 3))]
    (doseq [{:keys [expr]} top3]
      (println expr))))
