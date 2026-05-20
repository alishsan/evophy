
(ns evophy.core
  (:gen-class)
  (:require [clojure.edn :as edn]
            [clojure.walk :as walk]
            [emmy.env :as e]
            [taoensso.timbre :as timbre]))

(defn- rewrite-div-in-expr [expr]
  (walk/postwalk #(if (= % 'e/div) 'e// %) expr))

(defn- simplify-expr
  "Small algebraic simplifier to remove dead GP code before compile/score."
  [expr]
  (letfn [(simp [x]
            (if (and (sequential? x) (seq x))
              (let [[op & args] (map simp x)]
                (case op
                  +
                  (let [[a b] args]
                    (cond
                      (= a 0) b
                      (= b 0) a
                      (and (number? a) (number? b)) (+ a b)
                      :else (list '+ a b)))
                  -
                  (let [[a b] args]
                    (cond
                      (= a b) 0
                      (= b 0) a
                      (and (number? a) (number? b)) (- a b)
                      :else (list '- a b)))
                  *
                  (let [[a b] args]
                    (cond
                      (or (= a 0) (= b 0)) 0
                      (= a 1) b
                      (= b 1) a
                      (and (number? a) (number? b)) (* a b)
                      :else (list '* a b)))
                  e/square
                  (let [[a] args]
                    (if (number? a) (* a a) (list 'e/square a)))
                  e/sin
                  (let [[a] args]
                    (if (number? a) (Math/sin (double a)) (list 'e/sin a)))
                  e/cos
                  (let [[a] args]
                    (if (number? a) (Math/cos (double a)) (list 'e/cos a)))
                  e/sqrt
                  (let [[a] args]
                    (if (and (number? a) (not (neg? a)))
                      (Math/sqrt (double a))
                      (list 'e/sqrt a)))
                  e/div
                  (let [[a b] args]
                    (if (and (number? a) (number? b) (not (zero? b)))
                      (/ a b)
                      (list 'e/div a b)))
                  (cons op args)))
              x))]
    (simp expr)))

(defn grav2d-energy
  "H = (px² + py²)/(2m) − α/r,  r = √(qx² + qy²)."
  [m alpha {:keys [qx qy px py]}]
  (let [r (max (Math/sqrt (+ (* qx qx) (* qy qy))) 1e-12)]
    (+ (/ (+ (* px px) (* py py)) (* 2.0 m))
       (- (/ alpha r)))))

(defn- grav2d-deriv
  "Hamilton's equations: q̇ = ∂H/∂p, ṗ = −∂H/∂q."
  [m alpha {:keys [qx qy px py]}]
  (let [r2 (+ (* qx qx) (* qy qy))
        r3 (max (* r2 (Math/sqrt r2)) 1e-12)
        scale (/ (- alpha) r3)]
    {:dqx (/ px m)
     :dqy (/ py m)
     :dpx (* scale qx)
     :dpy (* scale qy)}))

(defn- symplectic-step [m alpha state dt]
  (let [half (* 0.5 dt)
        d1 (grav2d-deriv m alpha state)
        px' (+ (:px state) (* (:dpx d1) half))
        py' (+ (:py state) (* (:dpy d1) half))
        qx' (+ (:qx state) (* (:dqx d1) half))
        qy' (+ (:qy state) (* (:dqy d1) half))
        mid {:qx qx' :qy qy' :px px' :py py'}
        d2 (grav2d-deriv m alpha mid)
        qx'' (+ qx' (* (:dqx d2) dt))
        qy'' (+ qy' (* (:dqy d2) dt))
        late (assoc mid :qx qx'' :qy qy'')
        d3 (grav2d-deriv m alpha late)]
    {:qx qx'' :qy qy''
     :px (+ px' (* (:dpx d3) half))
     :py (+ py' (* (:dpy d3) half))}))

(defn generate-data [m alpha q0x q0y p0x p0y dt steps]
  (let [ic {:qx (double q0x) :qy (double q0y) :px (double p0x) :py (double p0y)}
        states (loop [s ic, n 0, acc []]
                 (if (> n steps)
                   acc
                   (recur (symplectic-step m alpha s dt) (inc n) (conj acc s))))]
    (mapv (fn [n s]
            (assoc s :t (* (double dt) n)
                   :energy (grav2d-energy m alpha s)))
          (range (count states))
          states)))

(defn scenario-data
  "Scenario map plus integrated :data trajectory (includes :m :alpha for fitness)."
  [scenario]
  (let [{:keys [m alpha q0x q0y p0x p0y dt steps]} scenario]
    (assoc scenario :data (vec (generate-data m alpha q0x q0y p0x p0y dt steps)))))

(def default-scenarios
  "2D gravitational ICs (|q| > 0, away from collision singularity)."
  [{:id :x-orbit   :m 1.0 :alpha 1.0 :q0x 2.5 :q0y 0.0 :p0x 0.0 :p0y 0.35 :dt 0.05 :steps 100}
   {:id :y-orbit   :m 1.0 :alpha 1.0 :q0x 0.0 :q0y 2.5 :p0x 0.35 :p0y 0.0 :dt 0.05 :steps 100}
   {:id :diagonal  :m 1.0 :alpha 1.0 :q0x 2.0 :q0y 1.5 :p0x -0.2 :p0y 0.3 :dt 0.05 :steps 100}
   {:id :heavy-m   :m 2.0 :alpha 1.0 :q0x 2.2 :q0y 1.0 :p0x 0.1 :p0y 0.25 :dt 0.04 :steps 100}
   {:id :strong-g  :m 1.0 :alpha 2.0 :q0x 2.8 :q0y 0.5 :p0x 0.05 :p0y 0.4 :dt 0.04 :steps 80}])

(defn scenarios->datasets
  [scenarios]
  (mapv scenario-data scenarios))

(def ^:private ic-vars '[q0x q0y p0x p0y])
(def ^:private state-vars '[qx qy px py])
(def ^:private param-vars '[m alpha])

(defn- compile-state-fn [expr]
  (binding [*ns* (the-ns 'evophy.core)]
    (eval `(fn [~'t ~'q0x ~'q0y ~'p0x ~'p0y ~'m ~'alpha]
             ~(-> expr simplify-expr rewrite-div-in-expr)))))

(defn- compile-rate-fn [expr]
  (binding [*ns* (the-ns 'evophy.core)]
    (eval `(fn [~'qx ~'qy ~'px ~'py ~'m ~'alpha]
             ~(-> expr simplify-expr rewrite-div-in-expr)))))

(defn- state-at-t [data t]
  (some (fn [s]
          (when (< (Math/abs (- (double (:t s)) (double t))) 1e-8)
            s))
        data))

(defn- all-horizon-times [data]
  (vec (distinct (map :t (filter #(pos? (double (:t %))) data)))))

(defn- initial-ics [data]
  (let [{:keys [qx qy px py]} (first data)]
    {:q0x (double qx) :q0y (double qy) :p0x (double px) :p0y (double py)}))

(defn- horizon-errors-analytical [fns {:keys [data m alpha]} horizon-times]
  (let [{:keys [q0x q0y p0x p0y]} (initial-ics data)
        md (double m)
        ad (double alpha)
        {qx-fn :qx qy-fn :qy px-fn :px py-fn :py} fns]
    (reduce
     (fn [acc t]
       (if-let [{:keys [qx qy px py]} (state-at-t data t)]
         (let [td (double t)
               pred-qx (double (qx-fn td q0x q0y p0x p0y md ad))
               pred-qy (double (qy-fn td q0x q0y p0x p0y md ad))
               pred-px (double (px-fn td q0x q0y p0x p0y md ad))
               pred-py (double (py-fn td q0x q0y p0x p0y md ad))]
           (-> acc
               (update :sq-q + (+ (e/square (- pred-qx qx))
                                  (e/square (- pred-qy qy))))
               (update :sq-p + (+ (e/square (- pred-px px))
                                  (e/square (- pred-py py))))
               (update :n inc)))
         acc))
     {:sq-q 0.0 :sq-p 0.0 :n 0}
     horizon-times)))

(defn- predict-step-rates
  "One Euler step from current (qx,qy,px,py) using evolved rate laws."
  [{:keys [dqx dqy dpx dpy]} {:keys [qx qy px py]} dt m alpha]
  (let [md (double m)
        ad (double alpha)
        dqx (double (dqx qx qy px py md ad))
        dqy (double (dqy qx qy px py md ad))
        dpx (double (dpx qx qy px py md ad))
        dpy (double (dpy qx qy px py md ad))
        dt (double dt)]
    {:qx (+ (double qx) (* dqx dt))
     :qy (+ (double qy) (* dqy dt))
     :px (+ (double px) (* dpx dt))
     :py (+ (double py) (* dpy dt))}))

(defn- one-step-errors-differential [rate-fns {:keys [data m alpha]}]
  (reduce
   (fn [acc [s0 s1]]
     (let [dt (- (double (:t s1)) (double (:t s0)))
           pred (predict-step-rates rate-fns s0 dt m alpha)]
       (-> acc
           (update :sq-q + (+ (e/square (- (:qx pred) (:qx s1)))
                              (e/square (- (:qy pred) (:qy s1)))))
           (update :sq-p + (+ (e/square (- (:px pred) (:px s1)))
                              (e/square (- (:py pred) (:py s1)))))
           (update :n inc))))
   {:sq-q 0.0 :sq-p 0.0 :n 0}
   (partition 2 1 data)))

(defn- rollout-errors-differential [rate-fns {:keys [data m alpha]}]
  (if (<= (count data) 1)
    {:sq-q 0.0 :sq-p 0.0 :n 0}
    (let [s0 (first data)]
      (loop [pred {:qx (:qx s0) :qy (:qy s0) :px (:px s0) :py (:py s0)}
             prev s0
             remaining (rest data)
             acc {:sq-q 0.0 :sq-p 0.0 :n 0}]
        (if-let [actual (first remaining)]
          (let [dt (- (double (:t actual)) (double (:t prev)))
                pred-next (predict-step-rates rate-fns pred dt m alpha)]
            (recur pred-next
                   actual
                   (rest remaining)
                   (-> acc
                       (update :sq-q + (+ (e/square (- (:qx pred-next) (:qx actual)))
                                          (e/square (- (:qy pred-next) (:qy actual)))))
                       (update :sq-p + (+ (e/square (- (:px pred-next) (:px actual)))
                                          (e/square (- (:py pred-next) (:py actual)))))
                       (update :n inc))))
          acc)))))

(defn expr-uses-symbol? [expr sym]
  (cond
    (= expr sym) true
    (coll? expr) (boolean (some #(expr-uses-symbol? % sym) expr))
    :else false))

(defn expr-uses-t? [expr]
  (expr-uses-symbol? expr 't))

(def analytical-expr-keys [:qx-expr :qy-expr :px-expr :py-expr])
(def differential-expr-keys [:dqx-expr :dqy-expr :dpx-expr :dpy-expr])

(def required-analytical-symbols (into ic-vars param-vars))
(def ^:private required-differential-symbols (into state-vars param-vars))

(defn- symbols-covered-across-exprs? [expr-keys ind required-syms]
  (every? (fn [sym]
            (some #(expr-uses-symbol? (get ind %) sym) expr-keys))
          required-syms))

(defn- inject-symbol-into-expr [expr sym]
  (list '+ expr sym))

(defn ensure-symbol-coverage
  "Add missing required symbols into random expr slots so genome-valid? can succeed."
  [ind expr-keys required-syms]
  (reduce
   (fn [ind sym]
     (if (symbols-covered-across-exprs? expr-keys ind [sym])
       ind
       (let [k (rand-nth expr-keys)]
         (update ind k #(inject-symbol-into-expr % sym)))))
   ind
   required-syms))

(defn ensure-analytical-uses-t [ind]
  (into ind
        (map (fn [k]
               (let [expr (get ind k)]
                 [k (if (expr-uses-t? expr) expr (list '+ 't expr))]))
             analytical-expr-keys)))

(defn- analytical-genome-valid? [ind]
  (and (every? #(expr-uses-t? (get ind %)) analytical-expr-keys)
       (symbols-covered-across-exprs? analytical-expr-keys ind required-analytical-symbols)))

(defn genome-valid?
  "Analytical: each expr uses t; ICs and (m, α) appear across trajectory laws.
   Differential: state coords and (m, α) appear across rate laws."
  [{:keys [strategy] :as ind}]
  (case strategy
    :analytical (analytical-genome-valid? ind)
    :differential (symbols-covered-across-exprs? differential-expr-keys ind
                                                required-differential-symbols)
    false))

(defn individual-genome-key [ind]
  (case (:strategy ind)
    :analytical (into [:analytical] (mapv ind analytical-expr-keys))
    :differential (into [:differential] (mapv ind differential-expr-keys))
    [(:strategy ind)]))

(defn- quant-double ^double [^double x]
  (double (/ (Math/round (* x 100000.0)) 100000.0)))

(defn- sample-trajectory-indices
  "Sparse indices along a trajectory for probe reuse (cheap, fixed across runs for same step counts)."
  [n]
  (if (zero? n)
    []
    (vec (distinct (filter #(< % n)
                            (if (= n 1)
                              [0]
                              [0 1 (quot n 2) (dec n)]))))))

(defn- differential-probes-from-dataset [dataset]
  (let [{:keys [data m alpha]} dataset
        n (count data)]
    (for [i (sample-trajectory-indices n)
          :let [s (nth data i)]]
      {:qx (double (:qx s))
       :qy (double (:qy s))
       :px (double (:px s))
       :py (double (:py s))
       :m (double m)
       :alpha (double alpha)})))

(defn- analytical-probes-from-dataset [dataset]
  (let [{:keys [data m alpha]} dataset
        {:keys [q0x q0y p0x p0y]} (initial-ics data)
        n (count data)]
    (for [i (sample-trajectory-indices n)
          :let [s (nth data i)]]
      {:t (double (:t s))
       :q0x (double q0x)
       :q0y (double q0y)
       :p0x (double p0x)
       :p0y (double p0y)
       :m (double m)
       :alpha (double alpha)})))

(defn build-behavior-probes
  "Probe tuples sampled from integrated scenarios — used to hash individuals by outputs, not syntax."
  [datasets]
  {:differential (vec (distinct (mapcat differential-probes-from-dataset datasets)))
   :analytical (vec (distinct (mapcat analytical-probes-from-dataset datasets)))})

(defn individual-behavior-key
  "Stable key from rounded model outputs on [[build-behavior-probes]]; nil if invalid or eval fails — use genome key then."
  [ind {:keys [differential analytical]}]
  (try
    (when (genome-valid? ind)
      (case (:strategy ind)
        :differential
        (when (seq differential)
          (let [dqx (compile-rate-fn (:dqx-expr ind))
                dqy (compile-rate-fn (:dqy-expr ind))
                dpx (compile-rate-fn (:dpx-expr ind))
                dpy (compile-rate-fn (:dpy-expr ind))]
            [:differential
             (vec
              (for [{:keys [qx qy px py m alpha]} differential]
                [(quant-double (double (dqx qx qy px py m alpha)))
                 (quant-double (double (dqy qx qy px py m alpha)))
                 (quant-double (double (dpx qx qy px py m alpha)))
                 (quant-double (double (dpy qx qy px py m alpha)))]))]))
        :analytical
        (when (seq analytical)
          (let [qx-fn (compile-state-fn (:qx-expr ind))
                qy-fn (compile-state-fn (:qy-expr ind))
                px-fn (compile-state-fn (:px-expr ind))
                py-fn (compile-state-fn (:py-expr ind))]
            [:analytical
             (vec
              (for [{:keys [t q0x q0y p0x p0y m alpha]} analytical]
                [(quant-double (double (qx-fn t q0x q0y p0x p0y m alpha)))
                 (quant-double (double (qy-fn t q0x q0y p0x p0y m alpha)))
                 (quant-double (double (px-fn t q0x q0y p0x p0y m alpha)))
                 (quant-double (double (py-fn t q0x q0y p0x p0y m alpha)))]))]))
        nil))
    (catch Exception _ nil)))

(defn take-distinct-by-behavior
  "Keep top n ranked individuals with unique [[individual-behavior-key]] (fallback: genome key)."
  [n ranked probes]
  (loop [seen #{} out [] xs (seq ranked)]
    (if (or (= (count out) n) (nil? xs))
      out
      (let [ind (first xs)
            k (or (individual-behavior-key ind probes)
                  (individual-genome-key ind))]
        (if (contains? seen k)
          (recur seen out (rest xs))
          (recur (conj seen k) (conj out ind) (rest xs)))))))

(defn- evaluate-analytical-predictions [ind dataset]
  (let [fns {:qx (compile-state-fn (:qx-expr ind))
             :qy (compile-state-fn (:qy-expr ind))
             :px (compile-state-fn (:px-expr ind))
             :py (compile-state-fn (:py-expr ind))}
        {:keys [data]} dataset
        times (all-horizon-times data)
        {:keys [sq-q sq-p n]} (horizon-errors-analytical fns dataset times)]
    (if (zero? n)
      {:strategy :analytical :n-horizons 0
       :mse-q Double/POSITIVE_INFINITY :mse-p Double/POSITIVE_INFINITY
       :mse Double/POSITIVE_INFINITY}
      (let [n (double n)]
        {:strategy :analytical :n-horizons (long n)
         :mse-q (/ sq-q n) :mse-p (/ sq-p n) :mse (/ (+ sq-q sq-p) n)}))))

(defn- evaluate-differential-predictions [ind dataset]
  (let [fns {:dqx (compile-rate-fn (:dqx-expr ind))
             :dqy (compile-rate-fn (:dqy-expr ind))
             :dpx (compile-rate-fn (:dpx-expr ind))
             :dpy (compile-rate-fn (:dpy-expr ind))}
        {:keys [sq-q sq-p n]} (one-step-errors-differential fns dataset)]
    (if (zero? n)
      {:strategy :differential :n-horizons 0
       :mse-q Double/POSITIVE_INFINITY :mse-p Double/POSITIVE_INFINITY
       :mse Double/POSITIVE_INFINITY}
      (let [n (double n)]
        {:strategy :differential :n-horizons (long n)
         :mse-q (/ sq-q n) :mse-p (/ sq-p n) :mse (/ (+ sq-q sq-p) n)}))))

(defn evaluate-predictions
  "Horizon error for an individual genome vs integrated trajectory (all t > 0).
   dataset is a scenario map with :data, :m, :alpha."
  [individual dataset]
  (case (:strategy individual)
    :analytical (evaluate-analytical-predictions individual dataset)
    :differential (evaluate-differential-predictions individual dataset)
    {:mse Double/POSITIVE_INFINITY}))

(def ops '[+ - * e/square e/sin e/cos e/div e/sqrt])
(def analytical-vars (vec (concat '[t] ic-vars param-vars)))
(def differential-vars (vec (concat state-vars param-vars)))
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
   'e/square {:arity 1 :fn e/square}
   'e/sqrt  {:arity 1 :fn e/sqrt}})

(defn- random-atom [vars]
  (rand-nth (concat vars constants)))

(def unary-ops '#{e/square e/sin e/cos e/sqrt})

(defn random-expression
  ([depth] (random-expression depth analytical-vars))
  ([depth vars]
   (if (or (zero? depth) (< (rand) 0.3))
     (random-atom vars)
     (let [op (rand-nth ops)]
       (if (unary-ops op)
         (list op (random-expression (dec depth) vars))
         (list op
               (random-expression (dec depth) vars)
               (random-expression (dec depth) vars)))))))

(defn- random-valid-individual [genome-fn expr-keys required-syms valid?-fn]
  (loop [n 0]
    (let [ind (-> (genome-fn)
                  (ensure-symbol-coverage expr-keys required-syms))]
      (cond
        (valid?-fn ind) ind
        (> n 500) ind
        :else (recur (inc n))))))

(defn random-analytical-individual []
  (random-valid-individual
   #(-> (into {:strategy :analytical}
              (map (fn [k] [k (random-expression 4 analytical-vars)])
                   analytical-expr-keys))
        (ensure-symbol-coverage analytical-expr-keys required-analytical-symbols)
        ensure-analytical-uses-t)
   analytical-expr-keys
   required-analytical-symbols
   analytical-genome-valid?))

(defn random-differential-individual []
  (random-valid-individual
   #(into {:strategy :differential}
          (map (fn [k] [k (random-expression 4 differential-vars)])
               differential-expr-keys))
   differential-expr-keys
   required-differential-symbols
   #(genome-valid? %)))

(defn random-individual []
  (if (< (rand) 0.5)
    (random-analytical-individual)
    (random-differential-individual)))

(def terminals (vec (concat '[t] ic-vars param-vars [(fn [] (rand-nth [0.5 1.0 2.0]))])))

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

(def ^:private fitness-rmse-floor 0.01)

(defn- fitness-from-errors
  [{:keys [sq-q sq-p n]} & {:keys [t-factor]}]
  (if (zero? n)
    0.0
    (let [n (double n)
          sum-sq (+ sq-q sq-p)]
      (if (Double/isNaN sum-sq)
        0.0
        (let [rmse (Math/sqrt (/ sum-sq n))]
          (if (Double/isNaN rmse)
            0.0
            (let [quality (/ 1.0 (+ rmse fitness-rmse-floor))]
              (* quality (or t-factor 1.0)))))))))

(defn calculate-analytical-fitness [ind dataset]
  (try
    (when (analytical-genome-valid? ind)
      (let [{:keys [data]} dataset
            fns {:qx (compile-state-fn (:qx-expr ind))
                 :qy (compile-state-fn (:qy-expr ind))
                 :px (compile-state-fn (:px-expr ind))
                 :py (compile-state-fn (:py-expr ind))}
            errors (horizon-errors-analytical fns dataset (all-horizon-times data))]
        (fitness-from-errors errors :t-factor 1.15)))
    (catch Exception _ 0)))

(defn calculate-differential-fitness [ind dataset]
  (try
    (when (genome-valid? ind)
      (let [fns {:dqx (compile-rate-fn (:dqx-expr ind))
                 :dqy (compile-rate-fn (:dqy-expr ind))
                 :dpx (compile-rate-fn (:dpx-expr ind))
                 :dpy (compile-rate-fn (:dpy-expr ind))}
            one-step (one-step-errors-differential fns dataset)
            rollout (rollout-errors-differential fns dataset)
            errors {:sq-q (+ (:sq-q one-step) (:sq-q rollout))
                    :sq-p (+ (:sq-p one-step) (:sq-p rollout))
                    :n (+ (:n one-step) (:n rollout))}]
        (fitness-from-errors errors)))
    (catch Exception _ 0)))

(defn calculate-fitness
  "dataset is a scenario map with :data, :m, :alpha."
  [individual dataset]
  (or (case (:strategy individual)
        :analytical (calculate-analytical-fitness individual dataset)
        :differential (calculate-differential-fitness individual dataset)
        nil)
      0))

(defn calculate-fitness-scenarios
  "Minimum per-scenario fitness—a law must fit every trajectory, not just on average."
  [individual datasets]
  (let [fits (mapv #(calculate-fitness individual %) datasets)]
    (if (empty? fits)
      0
      (apply min fits))))

(defn- normalize-expr
  "GP trees must be plain lists for eval/print; mutate used to leave LazySeq subtrees."
  [expr]
  (cond
    (instance? clojure.lang.LazySeq expr) (normalize-expr (doall expr))
    (and (sequential? expr) (not (map? expr)))
    (apply list (map normalize-expr expr))
    :else expr))

(defn mutate
  ([expr] (mutate expr analytical-vars))
  ([expr vars]
   (normalize-expr
    (if (< (rand) 0.2)
      (random-expression 2 vars)
      (if (coll? expr)
        (let [op (first expr)]
          (cons op (mapv #(if (coll? %) (mutate % vars) %) (rest expr))))
        expr)))))

(defn mutate-individual [ind]
  (if (< (rand) 0.2)
    (random-individual)
    (case (:strategy ind)
      :analytical (into {:strategy :analytical}
                        (map (fn [k] [k (mutate (get ind k) analytical-vars)])
                             analytical-expr-keys))
      :differential (into {:strategy :differential}
                          (map (fn [k] [k (mutate (get ind k) differential-vars)])
                               differential-expr-keys))
      (random-individual))))


(def default-population-file "data/population.edn")
(def checkpoint-version 7)

(defn- individual-for-save [ind]
  (let [norm (fn [k] [k (normalize-expr (get ind k))])]
    (case (:strategy ind)
      :analytical (into {:strategy :analytical} (map norm analytical-expr-keys))
      :differential (into {:strategy :differential} (map norm differential-expr-keys))
      (into {:strategy (:strategy ind)}
            (map norm (concat analytical-expr-keys differential-expr-keys))))))

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

(def default-mcts-simulations 64)
(def default-mcts-inject 5)

(defn parse-args
  [args]
  (loop [opts {:fresh? false
               :path default-population-file
               :generations 50
               :population-size 50
               :mcts-simulations default-mcts-simulations
               :mcts-inject default-mcts-inject
               :mcts? true}
         xs args]
    (if (empty? xs)
      opts
      (let [[a & more] xs]
        (case a
          "--fresh" (recur (assoc opts :fresh? true) more)
          "--no-mcts" (recur (assoc opts :mcts? false) more)
          "--mcts-simulations" (recur (assoc opts :mcts-simulations (Long/parseLong (first more))) (rest more))
          "--mcts-inject" (recur (assoc opts :mcts-inject (Long/parseLong (first more))) (rest more))
          "--population" (recur (assoc opts :path (first more)) (rest more))
          "--generations" (recur (assoc opts :generations (Long/parseLong (first more))) (rest more))
          "--population-size" (recur (assoc opts :population-size (Long/parseLong (first more))) (rest more))
          (throw (ex-info "Unknown argument (try --fresh --no-mcts --mcts-simulations N --mcts-inject N --population PATH --generations N --population-size N)"
                          {:arg a})))))))

(defn evolve-generation
  [population datasets population-size & {:keys [mcts? mcts-simulations mcts-inject]
                                          :or {mcts? true
                                               mcts-simulations default-mcts-simulations
                                               mcts-inject default-mcts-inject}}]
  (let [injected (when (and mcts? (pos? mcts-inject))
                   (require 'evophy.mcts)
                   ((ns-resolve 'evophy.mcts 'generation-inject)
                    datasets mcts-simulations mcts-inject))
        population (into population (or injected []))
        elite-n (max 1 (quot population-size 5))
        branch-factor (long (Math/ceil (/ population-size (double elite-n))))]
    (->> population
         (map (fn [ind]
                (assoc ind :fitness (calculate-fitness-scenarios ind datasets))))
         (sort-by :fitness #(compare %2 %1))
         (take elite-n)
         (mapcat (fn [parent]
                   (cons parent
                         (take (dec branch-factor)
                               (repeatedly #(mutate-individual parent))))))
         (take population-size)
         (vec))))

(defn -main [& args]
  (timbre/merge-config! {:min-level :warn})
  (let [{:keys [fresh? path generations population-size mcts? mcts-simulations mcts-inject]}
        (parse-args args)
        scenarios default-scenarios
        datasets (scenarios->datasets scenarios)
        {:keys [population generations-run resumed?]}
        (resolve-initial-population {:fresh? fresh?
                                     :path path
                                     :population-size population-size})
        initial (normalize-population-size population population-size)
        evolve-opts {:mcts? mcts? :mcts-simulations mcts-simulations :mcts-inject mcts-inject}
        final-pop (reduce (fn [pop _] (apply evolve-generation pop datasets population-size evolve-opts))
                          initial
                          (range generations))
        total-generations (+ generations-run generations)
        ranked (->> final-pop
                    (map (fn [ind]
                           (assoc ind :fitness (calculate-fitness-scenarios ind datasets))))
                    (sort-by :fitness #(compare %2 %1))
                    vec)
        behavior-probes (build-behavior-probes datasets)
        best (first ranked)
        top5 (take-distinct-by-behavior 5 ranked behavior-probes)]
    (save-population! path final-pop
                      :generations-run total-generations
                      :population-size population-size)
    (println (if resumed? "resumed from" "started new; saved to") path)
    (println "total generations (cumulative):" total-generations)
    (println "scenarios:" (count scenarios))
    (when mcts?
      (println "hybrid GP+MCTS:" mcts-inject "injections x" mcts-simulations "simulations/gen"))
    (when-not mcts?
      (println "hybrid: MCTS disabled (--no-mcts)"))
    (println "\nTop 5:")
    (doseq [[i {:keys [strategy fitness] :as ind}] (map-indexed vector top5)]
      (println (str "  #" (inc i) " " (name strategy) " fitness=" (format "%.4f" fitness)))
      (case strategy
        :analytical
        (doseq [k analytical-expr-keys]
          (println (str "    " (name k) ":" (pr-str (normalize-expr (get ind k))))))
        :differential
        (doseq [k differential-expr-keys]
          (println (str "    " (name k) ":" (pr-str (normalize-expr (get ind k))))))
        (println "    (unknown strategy)")))
    (println "\nBest scenario MSE:")
    (doseq [scenario scenarios]
      (let [dataset (scenario-data scenario)
            metrics (evaluate-predictions best dataset)]
        (println (str "  " (name (:id scenario)) " " (name (:strategy metrics))
                      " mse=" (format "%.6f" (:mse metrics))))))))
