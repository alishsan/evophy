
(ns evophy.core
  (:gen-class)
  (:require [clojure.edn :as edn]
            [clojure.walk :as walk]
            [emmy.env :as e]
            [evophy.analytical-blocks :as blocks]
            [evophy.coords :as coords]
            [evophy.primitives :as prims]
            [taoensso.timbre :as timbre]))

(defn- rewrite-div-in-expr [expr]
  (walk/postwalk #(if (= % 'e/div) 'e// %) expr))

(defn- rewrite-real-ops-in-expr
  "Replace Emmy ops with real-only helpers so sqrt/sin of invalid domains yield NaN, not Complex."
  [expr]
  (walk/postwalk
   (fn [x]
     (if (and (sequential? x) (seq x))
       (let [[op & args] x]
         (case op
           e/sqrt   (list 'real-sqrt (first args))
           e/sin    (list 'real-sin (first args))
           e/cos    (list 'real-cos (first args))
           e/sinh   (list 'real-sinh (first args))
           e/cosh   (list 'real-cosh (first args))
           e/square (list 'real-square (first args))
           e/div    (list 'real-div (first args) (second args))
           e//      (list 'real-div (first args) (second args))
           e/if     (list 'if (first args)
                           (list 'safe-double (second args))
                           (list 'safe-double (nth args 2)))
           x))
       x))
   expr))

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
                  (if (= (count args) 1)
                    (let [[a] args]
                      (if (number? a) (- a) (list '- a)))
                    (let [[a b] args]
                      (cond
                        (= a b) 0
                        (= b 0) a
                        (and (number? a) (number? b)) (- a b)
                        :else (list '- a b))))
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
                  e/sinh
                  (let [[a] args]
                    (if (number? a) (Math/sinh (double a)) (list 'e/sinh a)))
                  e/cosh
                  (let [[a] args]
                    (if (number? a) (Math/cosh (double a)) (list 'e/cosh a)))
                  e/sqrt
                  (let [[a] args]
                    (if (and (number? a) (not (neg? a)))
                      (Math/sqrt (double a))
                      (list 'e/sqrt a)))
                  e/div
                  (let [[a b] args]
                    (cond
                      (and (number? a) (number? b) (not (zero? b))) (/ a b)
                      (= a b) 1              ; a/a → 1 (catches qy/qy, px/px, etc.)
                      (= b 1) a              ; a/1 → a
                      :else (list 'e/div a b)))
                  e/if
                  (let [[test then else] args]
                    (list 'e/if (simp test) (simp then) (simp else)))
                  (cons op args)))
              x))]
    (let [match-cancel-mul (fn [x]
                              (when (and (sequential? x) (= (first x) '*))
                                (let [[_ a rhs] x]
                                  (when (and (sequential? rhs) (= (first rhs) '-))
                                    (let [[_ num y] rhs]
                                      (when (and (sequential? num)
                                                 (= (first num) 'e/div)
                                                 (= (nth num 2) a))
                                        (list '- (nth num 1) (list '* a y))))))))
          s (simp expr)]
      (cond
        (number? s) s
        (= s 0) 0
        ;; (+ a (- k a)) → k
        (and (sequential? s) (= (first s) '+))
        (let [[_ a b] s]
          (cond
            (and (sequential? b) (= (first b) '-)
                 (= (nth b 1) a)) (nth b 2)
            (and (sequential? a) (= (first a) '-)
                 (= (nth a 1) b)) (nth a 2)
            :else s))
        :else (or (match-cancel-mul s) s)))))

(defn grav2d-energy
  "H = (px² + py²)/(2m) − α/r,  r = √(qx² + qy²)."
  [m alpha {:keys [qx qy px py]}]
  (let [r (max (Math/sqrt (+ (* qx qx) (* qy qy))) 1e-12)]
    (+ (/ (+ (* px px) (* py py)) (* 2.0 m))
       (- (/ alpha r)))))

(def analytical-domains
  "Explicit validity of an analytical law: bound (E<0), unbound (E>0), or :any (legacy)."
  #{:bound :unbound :any})

(defn analytical-domain? [x]
  (contains? analytical-domains (keyword x)))

(defn scenario-regime
  "Orbit class from total energy at t=0: :bound, :unbound, or :parabolic (E≈0)."
  [{:keys [m alpha q0x q0y p0x p0y data]}]
  (let [{:keys [qx qy px py]}
        (if (seq data)
          (let [{:keys [qx qy px py]} (first data)]
            {:qx qx :qy qy :px px :py py})
          {:qx q0x :qy q0y :px p0x :py p0y})
        E (grav2d-energy m alpha {:qx qx :qy qy :px px :py py})]
    (cond
      (neg? (double E)) :bound
      (pos? (double E)) :unbound
      :else :parabolic)))

(defn analytical-law-domain
  "Domain tag on an analytical law; missing/invalid → :any (all scenarios)."
  [law]
  (let [d (:domain law)]
    (if (analytical-domain? d) (keyword d) :any)))

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

;; ── Exact Keplerian solver ────────────────────────────────────────────────
;; Unified complex oracle in evophy.primitives (M = u − e·sin u, n = (−2E)^(3/2)/(α√m)).

(defn- kepler-orbit-nan []
  (prims/kepler-orbit-nan))

(defn- kepler-orbit-at-t
  [m alpha q0x q0y p0x p0y t]
  (prims/kepler-orbit-at-t m alpha q0x q0y p0x p0y t))

(defn- kepler-state-at-t
  "Exact Keplerian (qx, qy, px, py) at time t via orbital elements.
   Handles elliptic (E < 0) and hyperbolic (E > 0).
   Returns nil for degenerate orbits (radial: |L| < 1e-8; near-parabolic: |e−1| < 1e-5)
   so callers can fall back to numerical integration."
  [m alpha q0x q0y p0x p0y t]
  (when-let [{:keys [qx qy px py]} (kepler-orbit-at-t m alpha q0x q0y p0x p0y t)]
    (when (and (Double/isFinite (double qx)) (Double/isFinite (double qy)))
      {:qx qx :qy qy :px px :py py})))

(defn generate-data [m alpha q0x q0y p0x p0y dt steps]
  (let [ic           {:qx (double q0x) :qy (double q0y) :px (double p0x) :py (double p0y)}
        use-kepler?  (some? (kepler-state-at-t m alpha q0x q0y p0x p0y 0.0))]
    (if use-kepler?
      ;; Exact analytic Keplerian states — no drift across any number of periods.
      (mapv (fn [n]
              (let [ti (* (double dt) n)
                    s  (kepler-state-at-t m alpha q0x q0y p0x p0y ti)]
                (assoc s :t ti :energy (grav2d-energy m alpha s))))
            (range (inc steps)))
      ;; Fallback: symplectic ODE integration for degenerate cases (radial / near-parabolic).
      (let [states (loop [s ic n 0 acc []]
                     (if (> n steps)
                       acc
                       (recur (symplectic-step m alpha s dt) (inc n) (conj acc s))))]
        (mapv (fn [n s]
                (assoc s :t (* (double dt) n)
                       :energy (grav2d-energy m alpha s)))
              (range (count states))
              states)))))

(defn scenario-data
  "Scenario map plus integrated :data trajectory (includes :m :alpha for fitness).
   Each point also carries polar (r, θ, p_r, p_θ); :chart is the preferred coordinate chart."
  [scenario]
  (let [{:keys [m alpha q0x q0y p0x p0y dt steps]} scenario
        base  (assoc scenario :data (vec (generate-data m alpha q0x q0y p0x p0y dt steps)))
        regime (scenario-regime base)]
    (-> base
        (coords/enrich-dataset)
        (assoc :regime regime))))

(defn- orbital-period
  "Keplerian period T = 2π √(m/α) · a^(3/2) for a bound orbit (E < 0).
   Returns nil for unbound (E ≥ 0)."
  [m alpha energy]
  (when (neg? (double energy))
    (let [a (/ (double alpha) (* -2.0 (double energy)))]
      (* 2.0 Math/PI (Math/sqrt (* (/ (double m) (double alpha)) a a a))))))

(def default-scenarios
  "Reference scenarios covering circular, low/high-eccentricity elliptic (axis-aligned
   and tilted), heavy-mass, strong-gravity, and hyperbolic orbits.
   Steps are auto-sized so each bound orbit is integrated for ≥ 1.2 full periods."
  (let [make (fn [sc]
               (let [{:keys [m alpha q0x q0y p0x p0y dt]} sc
                     e0    (grav2d-energy m alpha {:qx q0x :qy q0y :px p0x :py p0y})
                     T     (orbital-period m alpha e0)
                     steps (if T (int (Math/ceil (/ (* 1.2 T) (double dt)))) 200)]
                 (assoc sc :steps steps)))]
    ;; v_c(r) = √(α/mr).  Circular ICs: p0y = v_c, p0x = 0, q0y = 0.
    [(make {:id :circle         :m 1.0 :alpha 1.0 :q0x  2.0 :q0y 0.0 :p0x  0.0   :p0y 0.7071 :dt 0.10})
     ;; Low eccentricity (e ≈ 0.37) — moderate deviation from circular.
     (make {:id :ellipse-lo-e   :m 1.0 :alpha 1.0 :q0x  2.5 :q0y 0.0 :p0x  0.0   :p0y 0.50   :dt 0.05})
     ;; High eccentricity (e ≈ 0.60) — large r variation, strong force at periapsis.
     (make {:id :ellipse-hi-e   :m 1.0 :alpha 1.0 :q0x  4.0 :q0y 0.0 :p0x  0.0   :p0y 0.316  :dt 0.10})
     ;; Tilted ellipse — force has significant components in both x and y throughout.
     (make {:id :ellipse-tilted :m 1.0 :alpha 1.0 :q0x  1.5 :q0y 1.5 :p0x -0.35  :p0y 0.35   :dt 0.05})
     ;; Heavy mass — same geometry, different inertia; tests m dependence.
     (make {:id :heavy-m        :m 2.0 :alpha 1.0 :q0x  2.0 :q0y 0.0 :p0x  0.0   :p0y 0.50   :dt 0.04})
     ;; Stronger gravity — tests α dependence independently of m.
     (make {:id :strong-g       :m 1.0 :alpha 2.0 :q0x  2.0 :q0y 0.0 :p0x  0.0   :p0y 0.90   :dt 0.04})
     ;; Hyperbolic (E > 0) flyby — open trajectory, tests unbound regime.
     {:id :hyperbola :m 1.0 :alpha 1.0 :q0x -5.0 :q0y 2.5 :p0x 0.8 :p0y 0.0 :dt 0.10 :steps 200}]))

(defn scenarios->datasets
  [scenarios]
  (mapv scenario-data scenarios))

(defonce ^:private default-ref-datasets
  (delay (scenarios->datasets default-scenarios)))

(defn- reference-datasets []
  @default-ref-datasets)

(def default-scenario-bounds
  "Box in (m, α, q₀, p₀) for random scenario sampling; |q₀| ≥ :r-min.
   Wider ranges than before to include hyperbolic-energy scenarios in training."
  {:m [0.75 2.5]
   :alpha [0.75 2.5]
   :q0x [-4.0 4.0]
   :q0y [-4.0 4.0]
   :p0x [-0.8 0.8]
   :p0y [-0.8 0.8]
   :dt [0.04 0.10]
   :steps [100 250]
   :r-min 1.0})

(defn- uniform-sample
  "Uniform double in [lo, hi]; uses rng when non-nil, else clojure.core/rand."
  [rng lo hi]
  (if rng
    (+ lo (* (- hi lo) (.nextDouble ^java.util.Random rng)))
    (+ lo (* (- hi lo) (rand)))))

(defn- uniform-int-sample
  [rng lo hi]
  (+ lo (if rng
          (.nextInt ^java.util.Random rng (inc (- hi lo)))
          (rand-int (inc (- hi lo))))))

(defn random-scenario
  "Sample one scenario map inside [[default-scenario-bounds]] (|q₀| ≥ r-min)."
  [& {:keys [bounds rng id]}]
  (let [bounds (or bounds default-scenario-bounds)
        {:keys [m alpha q0x q0y p0x p0y dt steps r-min]} bounds
        r-min (or r-min 1.0)]
    (loop [attempt 0]
      (let [q0x (uniform-sample rng (first q0x) (second q0x))
            q0y (uniform-sample rng (first q0y) (second q0y))
            r (Math/sqrt (+ (* q0x q0x) (* q0y q0y)))]
        (if (>= r r-min)
          {:id (or id (keyword (str "rand-" (System/nanoTime))))
           :m (uniform-sample rng (first m) (second m))
           :alpha (uniform-sample rng (first alpha) (second alpha))
           :q0x q0x
           :q0y q0y
           :p0x (uniform-sample rng (first p0x) (second p0x))
           :p0y (uniform-sample rng (first p0y) (second p0y))
           :dt (uniform-sample rng (first dt) (second dt))
           :steps (uniform-int-sample rng (first steps) (second steps))}
          (if (> attempt 200)
            (throw (ex-info "Could not sample valid |q0| >= r-min" {:bounds bounds :r-min r-min}))
            (recur (inc attempt))))))))

(defn sample-random-scenarios
  "Sample n scenarios; optional :seed for reproducible draws."
  [n & {:keys [bounds seed]}]
  (let [rng (when (some? seed) (java.util.Random. (long seed)))]
    (mapv (fn [i] (random-scenario :bounds bounds :rng rng :id (keyword (str "rand-" i))))
          (range n))))

(defn make-fitness-context
  "Build a fitness evaluation context.

  :mode — :random (default) or :fixed (uses [[default-scenarios]]); used when
  :evaluation is :data-driven.
  :evaluation — :data-driven (trajectory fit) or :de-driven (analytical: 20% + full-orbit
  trajectory fit on integrated orbits; conserved: invariance along integrated orbits only —
  conserved formulas are not DE solutions). Differential is redundant in :de-driven.
  :sample-count — random scenarios per evaluation batch (default 32).
  :phase-samples — random phase-space points per batch when :de-driven (48).
  :aggregate — :min (worst case) or :percentile.
  :percentile — fraction in [0,1] or whole percent e.g. 10 for 10th percentile (default 0.1).
  :seed — optional; with :generation offset, stabilizes per-generation batches across a run."
  [& {:keys [mode evaluation sample-count phase-samples aggregate percentile seed]
      :or {mode :random
           evaluation :data-driven
           sample-count 32
           phase-samples 48
           aggregate :min
           percentile 0.1}}]
  (let [pct (if (<= percentile 1) (double percentile) (/ (double percentile) 100.0))]
    {:mode (keyword mode)
     :evaluation (keyword evaluation)
     :sample-count (long sample-count)
     :phase-samples (long phase-samples)
     :aggregate (keyword aggregate)
     :percentile pct
     :seed (when seed (long seed))
     :datasets (when (= :fixed (keyword mode))
                 (reference-datasets))}))

(defn datasets-for-fitness-context
  "Materialize scenario datasets for one evaluation (one GP generation when :generation is set)."
  [ctx & {:keys [generation]}]
  (case (:mode ctx)
    :fixed (:datasets ctx)
    :random
    (let [n (:sample-count ctx 32)
          seed (when (:seed ctx) (+ (long (:seed ctx)) (long (or generation 0))))
          scenarios (apply sample-random-scenarios n
                           (cond-> {:bounds default-scenario-bounds}
                             seed (assoc :seed seed)))]
      (scenarios->datasets scenarios))
    (throw (ex-info "Unknown fitness context mode" {:mode (:mode ctx)}))))

(def ^:private ic-vars '[q0x q0y p0x p0y])
(def ^:private state-vars '[qx qy px py])
(def ^:private param-vars '[m alpha])
;; Derived geometric variables available only to differential rate expressions.
;; r  = sqrt(qx²+qy²)  — orbital radius
;; r2 = qx²+qy²        — radius squared (avoids redundant sqrt when dividing by r²)
;; r3 = r³             — appears in gravitational force denominator
(def ^:private derived-vars '[r r2 r3])

(defn- safe-double
  "Coerce v to double. Returns NaN for any non-Number (e.g. Complex from sqrt of negative)."
  ^double [v]
  (if (instance? Number v) (double v) Double/NaN))

(defn- real-sqrt [a]
  (let [x (safe-double a)]
    (if (and (Double/isFinite x) (not (neg? x))) (Math/sqrt x) Double/NaN)))

(defn- real-sin [a] (Math/sin (safe-double a)))

(defn- real-cos [a] (Math/cos (safe-double a)))

(defn- real-sinh [a] (Math/sinh (safe-double a)))

(defn- real-cosh [a] (Math/cosh (safe-double a)))

;; Graded primitive delegates — visible to eval in compile-state-fn.
(defn mean-anomaly [n t m0] (prims/mean-anomaly n t m0))
(defn anomaly-from-M [e energy M] (prims/anomaly-from-M e energy M))
(defn periapsis-xp [semi-a ecc energy anomaly]
  (prims/periapsis-xp semi-a ecc energy anomaly))
(defn periapsis-yp [semi-b bh energy anomaly]
  (prims/periapsis-yp semi-b bh energy anomaly))
(defn lab-qx [cos-om sin-om xp yp] (prims/lab-qx cos-om sin-om xp yp))
(defn lab-qy [sin-om cos-om xp yp] (prims/lab-qy sin-om cos-om xp yp))

(defn- real-square [a]
  (let [x (safe-double a)] (* x x)))

(defn- real-div [a b]
  (let [x (safe-double a)
        y (safe-double b)]
    (if (and (Double/isFinite x) (Double/isFinite y) (not (zero? y)))
      (/ x y)
      Double/NaN)))

(defn- format-mse [mse]
  (cond
    (nil? mse) "n/a"
    :else
    (let [x (safe-double mse)]
      (cond
        (Double/isNaN x) "n/a"
        (Double/isFinite x) (format "%.6f" x)
        :else "Infinity"))))

(def ^:private max-expr-tree-size 140)
(def ^:private max-simplify-tree-size 80)

(defn- expr-tree-size-fast
  "Node count without normalize-expr; aborts early on pathological depth."
  ([expr] (expr-tree-size-fast expr 0))
  ([expr depth]
   (cond
     (> depth 48) (inc max-expr-tree-size)
     (not (coll? expr)) 1
     :else (+ 1 (reduce + 0 (map #(expr-tree-size-fast % (inc depth)) expr))))))

(defn- expr-too-large? [expr]
  (>= (expr-tree-size-fast expr) max-expr-tree-size))

(defn- simplify-if-small [expr]
  (if (>= (expr-tree-size-fast expr) max-simplify-tree-size)
    expr
    (simplify-expr expr)))

(defn- prepare-expr-for-compile [expr]
  (-> expr simplify-if-small rewrite-div-in-expr
      prims/rewrite-primitives-in-expr
      rewrite-real-ops-in-expr))

(def ^:private ^java.util.concurrent.ConcurrentHashMap compile-fn-cache
  (java.util.concurrent.ConcurrentHashMap.))

(defn- compile-nan-fn
  "Fallback when an expression is too large to compile safely."
  [& _] Double/NaN)

(defn- compile-cached
  "Memoize eval-compiled fns; oversized trees compile to a constant NaN."
  [kind expr body-fn]
  (if (expr-too-large? expr)
    compile-nan-fn
    (let [key (str kind \space (pr-str expr))]
      (if-some [hit (.get compile-fn-cache key)]
        hit
        (let [f (body-fn)]
          (.put compile-fn-cache key f)
          f)))))

(defn- compile-state-fn [expr]
  ;; Derived variables from initial conditions:
  ;;   r0/r02/r03  — initial orbital radius and its powers
  ;;   omega       — circular angular velocity sqrt(α/(m·r₀³)); exact for circular orbit
  ;;   omega-L     — true initial angular velocity L/(m·r₀²) = (q×p)/(m·r₀²); works for all orbits
  ;; Together these let the GP express Taylor corrections, circular rotation, and the
  ;; angular-momentum rotation approximation q(t) ≈ R(Ω_L·t)·q₀.
  (compile-cached :state expr
                  (fn []
                    (binding [*ns* (the-ns 'evophy.core)]
                      (eval `(fn [~'t ~'q0x ~'q0y ~'p0x ~'p0y ~'m ~'alpha]
                               (let [~'r02    (+ (* ~'q0x ~'q0x) (* ~'q0y ~'q0y))
                                     ~'r0     (Math/sqrt (max ~'r02 1e-24))
                                     ~'r03    (* ~'r0 ~'r02)
                                     ~'r05    (* ~'r03 ~'r02)
                                     ~'r06    (* ~'r03 ~'r03)
                                     ~'omega  (Math/sqrt (max 0.0 (/ ~'alpha (* ~'m ~'r03))))
                                     ~'omega-L (/ (- (* ~'q0x ~'p0y) (* ~'q0y ~'p0x))
                                                  (* ~'m ~'r02))
                                     ~'energy (grav2d-energy ~'m ~'alpha {:qx ~'q0x :qy ~'q0y :px ~'p0x :py ~'p0y})
                                     ~'kepler (kepler-orbit-at-t ~'m ~'alpha ~'q0x ~'q0y ~'p0x ~'p0y ~'t)
                                     ~'kepler-M0 (prims/M0-at-ic ~'m ~'alpha ~'q0x ~'q0y ~'p0x ~'p0y)
                                     ~'ecc       (:ecc ~'kepler)
                                     ~'cos-om    (:cos-om ~'kepler)
                                     ~'sin-om    (:sin-om ~'kepler)
                                     ~'semi-a    (:semi-a ~'kepler)
                                     ~'semi-b    (:semi-b ~'kepler)
                                     ~'bh        (:bh ~'kepler)
                                     ~'n-mean    (:n-mean ~'kepler)
                                     ~'kepler-M  (+ ~'kepler-M0 (* ~'n-mean ~'t))
                                     ~'kepler-u  (:kepler-u ~'kepler)
                                     ~'kepler-F  (:kepler-F ~'kepler)
                                     ~'kepler-xp (:kepler-xp ~'kepler)
                                     ~'kepler-yp (:kepler-yp ~'kepler)
                                     ~'kepler-vxp (:kepler-vxp ~'kepler)
                                     ~'kepler-vyp (:kepler-vyp ~'kepler)]
                                 (safe-double
                                  ~(prepare-expr-for-compile expr)))))))))

(defn- compile-rate-fn [expr]
  (compile-cached :rate expr
                  (fn []
                    (binding [*ns* (the-ns 'evophy.core)]
                      (eval `(fn [~'qx ~'qy ~'px ~'py ~'m ~'alpha]
                               (let [~'r2 (+ (* ~'qx ~'qx) (* ~'qy ~'qy))
                                     ~'r  (Math/sqrt (max ~'r2 1e-24))
                                     ~'r3 (* ~'r ~'r2)]
                                 (safe-double
                                  ~(prepare-expr-for-compile expr)))))))))

(defn- compile-conserved-fn
  "Compile conserved-quantity expression to a fn [qx qy px py m alpha] -> double.
   Derives r, r2, r3 from position — identical binding to compile-rate-fn."
  [expr]
  (compile-cached :conserved expr
                  (fn []
                    (binding [*ns* (the-ns 'evophy.core)]
                      (eval `(fn [~'qx ~'qy ~'px ~'py ~'m ~'alpha]
                               (let [~'r2 (+ (* ~'qx ~'qx) (* ~'qy ~'qy))
                                     ~'r  (Math/sqrt (max ~'r2 1e-24))
                                     ~'r3 (* ~'r ~'r2)]
                                 (safe-double
                                  ~(prepare-expr-for-compile expr)))))))))

(defn- compile-polar-state-fn [expr]
  "Analytical law in polar chart: (t, r₀, θ₀, p_r₀, p_θ₀, m, α) → scalar."
  (compile-cached :polar-state expr
                  (fn []
                    (binding [*ns* (the-ns 'evophy.core)]
                      (eval `(fn [~'t ~'r0 ~'theta0 ~'pr0 ~'ptheta0 ~'m ~'alpha]
                               (let [~'r02    (* ~'r0 ~'r0)
                                     ~'r03    (* ~'r02 ~'r0)
                                     ~'omega  (Math/sqrt (max 0.0 (/ ~'alpha (* ~'m ~'r03))))
                                     ~'omega-L (/ ~'ptheta0 (* ~'m ~'r02))
                                     ~'energy (+ (/ (+ (* ~'pr0 ~'pr0)
                                                      (/ (* ~'ptheta0 ~'ptheta0) ~'r02))
                                                  (* 2.0 ~'m))
                                                 (- (/ ~'alpha ~'r0)))]
                                 (safe-double
                                  ~(prepare-expr-for-compile expr)))))))))

(defn- compile-polar-conserved-fn [expr]
  "Conserved law in polar chart: C(r, θ, p_r, p_θ, m, α)."
  [expr]
  (compile-cached :polar-conserved expr
                  (fn []
                    (binding [*ns* (the-ns 'evophy.core)]
                      (eval `(fn [~'r ~'theta ~'pr ~'ptheta ~'m ~'alpha]
                               (let [~'r2 (* ~'r ~'r)
                                     ~'r3 (* ~'r ~'r2)]
                                 (safe-double
                                  ~(prepare-expr-for-compile expr)))))))))

(defn- compile-conserved-fn-for-chart
  [chart expr]
  (if (= :polar (keyword chart))
    (compile-polar-conserved-fn expr)
    (compile-conserved-fn expr)))

(defn- compile-analytical-fns [ind chart]
  (if (= :polar (keyword chart))
    {:r     (compile-polar-state-fn (:r-expr ind))
     :theta (compile-polar-state-fn (:theta-expr ind))
     :pr    (compile-polar-state-fn (:pr-expr ind))
     :ptheta (compile-polar-state-fn (:ptheta-expr ind))}
    {:qx (compile-state-fn (:qx-expr ind))
     :qy (compile-state-fn (:qy-expr ind))
     :px (compile-state-fn (:px-expr ind))
     :py (compile-state-fn (:py-expr ind))}))

(defn- state-at-t [data t]
  (some (fn [s]
          (when (< (Math/abs (- (double (:t s)) (double t))) 1e-8)
            s))
        data))

(defn- all-horizon-times [data]
  (vec (distinct (map :t (filter #(pos? (double (:t %))) data)))))

(def analytical-gp-horizon-frac
  "Minimum training horizon: first 20% of each reference orbit."
  0.20)

(def analytical-full-horizon-frac
  "Full-orbit training horizon (all integrated samples)."
  1.0)

(def analytical-training-horizon-fracs
  "GP / DE-driven analytical fitness uses min fitness across these horizons."
  [analytical-gp-horizon-frac analytical-full-horizon-frac])

(defn- horizon-fraction-times
  "Time samples covering the first `frac` of a trajectory; frac >= 1 uses all samples."
  ([data] (horizon-fraction-times data analytical-gp-horizon-frac))
  ([data frac]
   (let [ts (all-horizon-times data)]
     (if (>= (double frac) 1.0)
       ts
       (vec (take (max 1 (int (Math/ceil (* (double frac) (count ts))))) ts))))))

(def repair-horizon-frac
  "MCTS adversarial repair uses the same dual training horizons as GP (via one-arg fitness)."
  analytical-gp-horizon-frac)

(defn- initial-ics [data]
  (let [{:keys [qx qy px py]} (first data)]
    {:q0x (double qx) :q0y (double qy) :p0x (double px) :p0y (double py)}))

(defn- horizon-errors-polar-analytical [fns {:keys [data m alpha]} horizon-times]
  (let [{:keys [r0 theta0 pr0 ptheta0]} (coords/initial-polar-ics data)
        md (double m) ad (double alpha)
        {r-fn :r theta-fn :theta pr-fn :pr ptheta-fn :ptheta} fns]
    (reduce
     (fn [acc t]
       (if-let [{:keys [r theta pr ptheta]} (state-at-t data t)]
         (let [td (double t)
               pred-r (safe-double (r-fn td r0 theta0 pr0 ptheta0 md ad))
               pred-th (safe-double (theta-fn td r0 theta0 pr0 ptheta0 md ad))
               pred-pr (safe-double (pr-fn td r0 theta0 pr0 ptheta0 md ad))
               pred-pt (safe-double (ptheta-fn td r0 theta0 pr0 ptheta0 md ad))]
           (-> acc
               (update :sq-q + (+ (e/square (- pred-r r))
                                  (e/square (- pred-th theta))))
               (update :sq-p + (+ (e/square (- pred-pr pr))
                                  (e/square (- pred-pt ptheta))))
               (update :n inc)))
         acc))
     {:sq-q 0.0 :sq-p 0.0 :n 0}
     horizon-times)))

(defn- horizon-errors-analytical [fns dataset horizon-times]
  (if (= :polar (:chart dataset))
    (horizon-errors-polar-analytical fns dataset horizon-times)
    (let [{:keys [data m alpha]} dataset
          {:keys [q0x q0y p0x p0y]} (initial-ics data)
          md (double m)
          ad (double alpha)
          {qx-fn :qx qy-fn :qy px-fn :px py-fn :py} fns]
      (reduce
       (fn [acc t]
         (if-let [{:keys [qx qy px py]} (state-at-t data t)]
           (let [td (double t)
                 pred-qx (safe-double (qx-fn td q0x q0y p0x p0y md ad))
                 pred-qy (safe-double (qy-fn td q0x q0y p0x p0y md ad))
                 pred-px (safe-double (px-fn td q0x q0y p0x p0y md ad))
                 pred-py (safe-double (py-fn td q0x q0y p0x p0y md ad))]
             (-> acc
                 (update :sq-q + (+ (e/square (- pred-qx qx))
                                      (e/square (- pred-qy qy))))
                 (update :sq-p + (+ (e/square (- pred-px px))
                                      (e/square (- pred-py py))))
                 (update :n inc)))
           acc))
       {:sq-q 0.0 :sq-p 0.0 :n 0}
       horizon-times))))

(defn- symplectic-step-rates
  "Störmer-Verlet (symplectic) step using evolved rate laws.
   Treats dqx/dqy as velocity equations and dpx/dpy as force equations,
   mirroring the integrator used to generate reference data.
   This ensures that even the correct Newtonian ODE accumulates near-zero
   rollout error, rather than the O(N·dt) drift that plain Euler produces."
  [{:keys [dqx dqy dpx dpy]} {:keys [qx qy px py]} dt m alpha]
  (let [md  (double m)
        ad  (double alpha)
        dt  (double dt)
        h   (* 0.5 dt)
        ;; Half-step momentum update
        fpx (safe-double (dpx qx qy px py md ad))
        fpy (safe-double (dpy qx qy px py md ad))
        px1 (+ (double px) (* fpx h))
        py1 (+ (double py) (* fpy h))
        ;; Full-step position update using mid-step velocity
        vqx (safe-double (dqx qx qy px1 py1 md ad))
        vqy (safe-double (dqy qx qy px1 py1 md ad))
        qx2 (+ (double qx) (* vqx dt))
        qy2 (+ (double qy) (* vqy dt))
        ;; Second half-step momentum update at new position
        fpx2 (safe-double (dpx qx2 qy2 px1 py1 md ad))
        fpy2 (safe-double (dpy qx2 qy2 px1 py1 md ad))]
    {:qx qx2 :qy qy2
     :px (+ px1 (* fpx2 h))
     :py (+ py1 (* fpy2 h))}))

(defn- one-step-errors-differential [rate-fns {:keys [data m alpha]}]
  (reduce
   (fn [acc [s0 s1]]
     (let [dt (- (double (:t s1)) (double (:t s0)))
           pred (symplectic-step-rates rate-fns s0 dt m alpha)]
       (-> acc
           (update :sq-q + (+ (e/square (- (:qx pred) (:qx s1)))
                              (e/square (- (:qy pred) (:qy s1)))))
           (update :sq-p + (+ (e/square (- (:px pred) (:px s1)))
                              (e/square (- (:py pred) (:py s1)))))
           (update :n inc))))
   {:sq-q 0.0 :sq-p 0.0 :n 0}
   (partition 2 1 data)))

(defn- rollout-errors-differential [rate-fns {:keys [data m alpha]}]
  ;; Uses symplectic integration for rollout so that the correct Newtonian ODE
  ;; accumulates near-zero error, not O(N·dt) Euler drift.
  (if (<= (count data) 1)
    {:sq-q 0.0 :sq-p 0.0 :n 0}
    (let [s0 (first data)]
      (loop [pred      {:qx (:qx s0) :qy (:qy s0) :px (:px s0) :py (:py s0)}
             prev      s0
             remaining (rest data)
             acc       {:sq-q 0.0 :sq-p 0.0 :n 0}]
        (if-let [actual (first remaining)]
          (let [dt         (- (double (:t actual)) (double (:t prev)))
                pred-next  (symplectic-step-rates rate-fns pred dt m alpha)]
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

(def ^:private kepler-time-proxy-syms
  "Kepler derived vars computed from t in compile-state-fn — satisfy trajectory-in-t checks."
  '#{kepler-u kepler-F kepler-M kepler-xp kepler-yp kepler-vxp kepler-vyp})

(defn expr-uses-t? [expr]
  (or (expr-uses-symbol? expr 't)
      (some #(expr-uses-symbol? expr %) kepler-time-proxy-syms)))

(def analytical-expr-keys    [:qx-expr :qy-expr :px-expr :py-expr])
(def differential-expr-keys  [:dqx-expr :dqy-expr :dpx-expr :dpy-expr])
(def conserved-expr-key      :c-expr)

;; ── Multi-law individuals ─────────────────────────────────────────────────────
;; An individual may bundle several formulas, e.g. analytical q(t),p(t) plus a
;; conserved C(q,p).  They are independent roles:
;;   analytical — candidate solution of the motion equations (q,p as functions of t)
;;   conserved  — scalar invariant obeyed along physical motion (not a DE solution)
;; Legacy checkpoints with a single top-level :strategy are wrapped on load.

(defn- law-kind [law]
  (or (:kind law) (:strategy law)))

(defn- law-from-legacy [ind]
  (cond-> (-> ind (dissoc :strategy) (assoc :kind (:strategy ind)))
    (analytical-domain? (:domain ind))
    (assoc :domain (keyword (:domain ind)))))

(defn- law->legacy [law]
  (let [k (law-kind law)]
    (cond-> (-> law (assoc :strategy k) (dissoc :kind))
      (and (= k :analytical) (not= (analytical-law-domain law) :any))
      (assoc :domain (analytical-law-domain law)))))

(defn individual-laws
  "All formula groups in an individual (legacy single-strategy maps become one law)."
  [ind]
  (cond
    (seq (:laws ind)) (vec (:laws ind))
    (:strategy ind)     [(law-from-legacy ind)]
    :else               []))

(defn laws-individual [laws]
  {:laws (vec laws)})

(defn individual-law-kinds [ind]
  (mapv law-kind (individual-laws ind)))

(defn composite-individual? [ind]
  (> (count (individual-laws ind)) 1))

(defn primary-strategy-label [ind]
  (cond
    (composite-individual? ind) :composite
    :else (or (some law-kind (individual-laws ind)) :unknown)))

(defn first-law-legacy [ind]
  (law->legacy (if (seq (:laws ind)) (first (:laws ind)) ind)))

(defn- individual-with-law [ind law]
  "Preserve :laws wrapper when mutating/crossing composite individuals."
  (let [law (if (:kind law) law (law-from-legacy law))]
    (if (seq (:laws ind))
      (assoc ind :laws [law])
      (law->legacy law))))

(defn- strip-differential-laws [ind]
  (let [laws (vec (remove #(= :differential (law-kind %)) (individual-laws ind)))]
    (when (seq laws)
      (laws-individual laws))))

(defn- migrate-individual-to-laws [ind]
  (if (seq (:laws ind))
    ind
    (when (:strategy ind)
      (laws-individual [(law-from-legacy ind)]))))

;; Analytical expressions must reference the initial position (q0x, q0y) so they aren't
;; trivially constant.  We don't mandate p0x/p0y or alpha explicitly: the circular solution
;; is fully determined by position and omega = sqrt(alpha/(m*r03)), so alpha is implicit.
;; Requiring m ensures the GP uses mass (affects the circular frequency).
(def required-analytical-symbols   '[q0x q0y m])
(def ^:private required-differential-symbols (into state-vars param-vars))

(def ^:private ic-proxy-syms
  "Orbit elements derived from ICs in compile-state-fn — cover q0x/q0y mandate for Kepler laws."
  '#{ecc cos-om sin-om semi-a semi-b bh n-mean kepler-M0
    kepler-u kepler-F kepler-xp kepler-yp kepler-vxp kepler-vyp kepler-M
    r0 r02 r03 omega omega-L})

(defn- sym-covered-in-expr? [expr sym]
  (or (expr-uses-symbol? expr sym)
      (and (#{'q0x 'q0y} sym)
           (some #(expr-uses-symbol? expr %) ic-proxy-syms))
      (and (#{'qx 'qy} sym)
           (some #(expr-uses-symbol? expr %) derived-vars))))

(defn- symbols-covered-across-exprs? [expr-keys ind required-syms]
  (every? (fn [sym]
            (some (fn [k]
                    (sym-covered-in-expr? (get ind k) sym))
                  expr-keys))
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

(def ^:dynamic *preferred-law-chart*
  "Default chart for new analytical laws and validation when a law has no explicit :chart."
  :cartesian)

(defn- infer-law-chart [law]
  (or (coords/law-chart law)
      (when (coords/analytical-exprs-present? law :polar) :polar)
      (when (coords/analytical-exprs-present? law :cartesian) :cartesian)
      *preferred-law-chart*))

(defn- expr-energy-regime-test? [expr]
  (and (sequential? expr)
       (= 'neg? (first expr))
       (= 'energy (second expr))))

(defn- expr-regime-branch? [expr]
  (and (sequential? expr)
       (= 'e/if (first expr))
       (= 3 (count (rest expr)))
       (expr-energy-regime-test? (second expr))))

(declare normalize-expr)

(defn- expr-all-eifs-use-energy-test? [expr]
  "Every e/if in the tree must branch on (neg? energy) — rejects (e/if t …) and (e/if (neg? q0y) …)."
  (let [expr (normalize-expr expr)]
    (every? (fn [node]
              (or (not (sequential? node))
                  (not= 'e/if (first node))
                  (expr-energy-regime-test? (second node))))
            (tree-seq sequential? seq expr))))

(defn analytical-strict-energy-branches?
  "Both-regimes validity: every analytical slot is (e/if (neg? energy) … …) with no bogus e/if tests."
  [law]
  (let [leg   (if (:strategy law) law (law->legacy law))
        chart (infer-law-chart leg)
        keys  (coords/analytical-expr-keys-for-chart chart)]
    (boolean
     (and (seq keys)
          (every? (fn [k]
                    (when-let [ex (get leg k)]
                      (and (expr-regime-branch? ex)
                           (expr-all-eifs-use-energy-test? ex))))
                  keys)))))

(defn- regime-arm-expr
  "Bound => then-branch; unbound (and parabolic) => else-branch of (e/if (neg? energy) … …)."
  [expr regime]
  (let [ex (normalize-expr expr)]
    (if (expr-regime-branch? ex)
      (let [args (vec (rest ex))
            then-branch (nth args 1)
            else-branch (nth args 2)]
        (if (= regime :bound) then-branch else-branch))
      ex)))

(defn- dataset-regime-for-fitness [dataset]
  (let [r (or (:regime dataset) (scenario-regime dataset))]
    (if (= r :bound) :bound :unbound)))

(defn analytical-law-regime-slice
  "Law with each slot replaced by its bound or unbound arm (no outer e/if)."
  [law regime]
  (let [leg   (if (:strategy law) law (law->legacy law))
        chart (infer-law-chart leg)
        keys  (coords/analytical-expr-keys-for-chart chart)]
    (into leg
          (map (fn [k] [k (regime-arm-expr (get leg k) regime)])
               keys))))

(defn- analytical-slice-valid? [ind chart]
  (let [keys     (coords/analytical-expr-keys-for-chart chart)
        required (coords/required-analytical-symbols-for-chart chart)]
    (and (every? #(expr-uses-t? (get ind %)) keys)
         (symbols-covered-across-exprs? keys ind required))))

(defn analytical-branches-on-energy?
  "True when any analytical slot uses (e/if (neg? energy) bound-branch unbound-branch)."
  [law]
  (let [leg   (if (:strategy law) law (law->legacy law))
        chart (infer-law-chart leg)
        keys  (coords/analytical-expr-keys-for-chart chart)]
    (boolean (some #(when-let [ex (get leg %)] (expr-regime-branch? ex)) keys))))

(def ^:private cartesian-analytical-expr-keys
  [:qx-expr :qy-expr :px-expr :py-expr])

(def ^:private unbound-linear-analytical-exprs
  "First-order Taylor — default unbound arm when wrapping bound seeds for both-regimes."
  '{:qx-expr (+ q0x (* (e/div p0x m) t))
    :qy-expr (+ q0y (* (e/div p0y m) t))
    :px-expr (+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
    :py-expr (+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))})

(defn- wrap-analytical-energy-branches
  "Each slot becomes (e/if (neg? energy) bound-expr unbound-expr)."
  [ind & {:keys [unbound-exprs] :or {unbound-exprs unbound-linear-analytical-exprs}}]
  (assoc (into ind
               (map (fn [k]
                      [k (list 'e/if '(neg? energy) (get ind k) (get unbound-exprs k))])
                    cartesian-analytical-expr-keys))
         :domain :any))

(def ^:dynamic *both-regimes?*
  "When true, analytical fitness/MSE use every scenario; genomes must branch on (neg? energy)."
  false)

(def ^:dynamic *template-unbound-arms?*
  "When true with both-regimes: lock unbound arms to template; evolve bound arms via q/p pair mutations."
  false)

(def ^:dynamic *template-conic-unbound?*
  "When true with template-unbound + analytical-blocks: unbound template is :hyperbola-conic (not Taylor)."
  true)

(def ^:dynamic *analytical-blocks?*
  "When true, inject Emmy-validated analytical catalog blocks (circle, ellipse, Taylor, …)."
  false)

(def ^:dynamic *primitive-tier*
  "Graded pure-function extension: 0 = algebra only; 1..5 unlock Kepler pipeline stages."
  0)

(def primitive-use-rate
  "Fraction of random GP trees that try one graded primitive call."
  0.15)

(defn current-primitive-tier [] (long *primitive-tier*))

(defn unlocked-gp-primitive-ops
  "Primitive ops available at the current tier (for MCTS / mutation)."
  [tier]
  (prims/unlocked-primitive-ops tier))

(defn template-unbound-arms-active? []
  (and *template-unbound-arms?* *both-regimes?*))

(defn- analytical-blocks-active? []
  (and *analytical-blocks?* *both-regimes?*))

(defn template-conic-unbound-active? []
  (and *template-conic-unbound?* (template-unbound-arms-active?) (analytical-blocks-active?)))

(def ^:private bound-arm-pairs
  "Hamiltonian-coupled slot pairs for coordinated bound-arm mutation."
  [[:qx-expr :px-expr] [:qy-expr :py-expr]])

(defn unbound-arm-expr
  "Template unbound arm for one analytical slot (Taylor or :hyperbola-conic)."
  [expr-key]
  (if-let [conic (when (template-conic-unbound-active?)
                   (blocks/hyperbola-conic-slots))]
    (get conic expr-key (get unbound-linear-analytical-exprs expr-key))
    (get unbound-linear-analytical-exprs expr-key)))

(defn slot-with-bound-arm
  "Wrap bound-arm in (e/if (neg? energy) bound Taylor-unbound)."
  [expr-key bound-expr]
  (list 'e/if '(neg? energy) bound-expr (unbound-arm-expr expr-key)))

(defn template-unbound-arms-law
  "Re-lock every unbound arm to Taylor while preserving bound arms."
  [leg]
  (if (and *template-unbound-arms?* *both-regimes?* (= :analytical (:strategy leg)))
    (into leg
          (map (fn [k]
                 (let [bound (regime-arm-expr (get leg k) :bound)]
                   [k (slot-with-bound-arm k bound)]))
               cartesian-analytical-expr-keys))
    leg))

(defn bound-arm-expr
  "First-order Taylor bound arm — pairs with unbound catalog blocks."
  [expr-key]
  (get blocks/taylor-bound-slots expr-key))

(defn slot-with-unbound-arm
  "Wrap unbound-arm in (e/if (neg? energy) Taylor-bound unbound)."
  [expr-key unbound-expr]
  (list 'e/if '(neg? energy) (bound-arm-expr expr-key) unbound-expr))

(defn- wrap-kepler-conic-strict
  "Both-regimes strict shell: :ellipse-conic bound + :hyperbola-conic unbound."
  [entry]
  (into {:strategy :analytical :domain :any}
        (map (fn [k]
               [k (list 'e/if '(neg? energy)
                        (get (:bound-slots entry) k)
                        (get (:unbound-slots entry) k))])
             cartesian-analytical-expr-keys)))

(defn- wrap-catalog-entry-strict
  "Strict e/if shell for one catalog entry (bound → template unbound; unbound → Taylor bound)."
  [entry]
  (cond
    (= :kepler-conic (:id entry))
    (wrap-kepler-conic-strict entry)

    :else
    (let [regime (:regime entry :bound)
          slots  (:slots entry)
          leg    (into {:strategy :analytical
                        :domain (if *both-regimes?* :any (blocks/regime->domain regime))}
                       slots)]
      (cond
        (not *both-regimes?*)
        leg

        (= regime :unbound)
        (into leg
              (map (fn [k] [k (slot-with-unbound-arm k (get slots k))])
                   cartesian-analytical-expr-keys))

        (template-unbound-arms-active?)
        (template-unbound-arms-law
         (into leg
               (map (fn [k] [k (slot-with-bound-arm k (get slots k))])
                    cartesian-analytical-expr-keys)))

        :else
        (wrap-analytical-energy-branches leg)))))

(defn- catalog-law-as-individual []
  (if (and (template-conic-unbound-active?) (blocks/kepler-conic-valid?) (< (rand) 0.40))
    (wrap-kepler-conic-strict (blocks/kepler-conic-entry))
    (wrap-catalog-entry-strict (blocks/random-catalog-entry))))

(defn- strict-catalog-law-from-entry [entry]
  (wrap-catalog-entry-strict entry))

(defn- apply-catalog-unbound-pair-to-law
  "Re-wrap one (q,p) pair: preserve bound arm, lock unbound to hyperbola-conic template."
  [leg]
  (let [pair (rand-nth (blocks/unbound-pairs))]
    (into leg
          (map (fn [slot]
                 [slot (slot-with-bound-arm slot (regime-arm-expr (get leg slot) :bound))])
               pair))))

(defn- apply-catalog-pair-to-law [leg]
  (let [entry   (blocks/random-catalog-entry)
        pair    (rand-nth (if (= :unbound (:regime entry))
                          (blocks/unbound-pairs)
                          (blocks/bound-pairs)))
        grafted (blocks/graft-catalog-pair leg :entry entry :pair pair)]
    (if (= :unbound (:regime entry))
      (into grafted
            (map (fn [slot]
                   [slot (slot-with-unbound-arm slot (get grafted slot))])
                 pair))
      (if (template-unbound-arms-active?)
        (into grafted
              (map (fn [slot]
                     [slot (slot-with-bound-arm slot (get grafted slot))])
                   pair))
        grafted))))

(defn- per-regime-arm-fitness? [ind]
  (and *both-regimes?* (analytical-strict-energy-branches? ind)))

(defn domain-applicable?
  "True when an analytical law should be scored on this scenario/dataset.
   Laws with (e/if (neg? energy) … …) apply to every scenario; :domain tags restrict.
   With *both-regimes?*, every scenario applies (domain tags ignored)."
  [law scenario-or-dataset]
  (or *both-regimes?*
      (analytical-branches-on-energy? law)
      (let [ld (analytical-law-domain law)
            sr (or (:regime scenario-or-dataset)
                 (scenario-regime scenario-or-dataset))]
        (or (= ld :any)
            (= ld sr)
            (and (= ld :unbound) (= sr :parabolic))))))

(defn- datasets-for-analytical-law [law datasets]
  (filterv #(domain-applicable? law %) datasets))

(defn- analytical-genome-valid? [ind]
  (let [chart (infer-law-chart ind)
        keys  (coords/analytical-expr-keys-for-chart chart)
        required (coords/required-analytical-symbols-for-chart chart)
        tier  (long *primitive-tier*)]
    (and (every? #(expr-uses-t? (get ind %)) keys)
         (symbols-covered-across-exprs? keys ind required)
         (every? #(prims/expr-primitive-valid? (get ind %) tier) keys)
         (or (not *both-regimes?*)
             (analytical-strict-energy-branches? ind)))))

(declare conserved-vals-for-dataset conserved-dot-vals-for-dataset)

(defn- conserved-per-traj-means
  "Per-trajectory mean of conserved expression over datasets."
  [c-fn datasets]
  (mapv (fn [ds]
          (let [vals   (conserved-vals-for-dataset c-fn ds)
                finite (filterv #(Double/isFinite %) vals)
                n      (count finite)]
            (when (pos? n) (/ (reduce + finite) n))))
        datasets))

(defn- conserved-expr-chart
  "Infer coordinate chart from symbols used in C — polar when pr/pθ appear."
  [ind]
  (let [expr (get ind conserved-expr-key)
        sexpr (when (some? expr) (simplify-expr expr))]
    (if (and sexpr (some #(expr-uses-symbol? sexpr %) '[pr ptheta]))
      :polar
      :cartesian)))

(defn- conserved-chart-compatible? [law chart]
  (or (coords/law-chart law)
      (= chart (conserved-expr-chart law))))

(defn- conserved-mean-for-dataset [ind ds]
  (let [expr-chart (conserved-expr-chart ind)
        c-fn  (compile-conserved-fn-for-chart expr-chart (:c-expr ind))
        vals  (conserved-vals-for-dataset c-fn (assoc ds :chart expr-chart))
        finite (filterv #(Double/isFinite %) vals)
        n      (count finite)]
    (when (pos? n) (/ (reduce + finite) n))))

(defn- ic-group-variation-ok?
  "True when at least one (m, α) group has ≥2 scenarios whose per-trajectory
   means differ by >5% CoV — the expression must depend on initial conditions,
   not just parameters.  Returns false when no qualifying groups exist."
  [ind datasets]
  (let [expr-chart (conserved-expr-chart ind)
        datasets   (filterv #(= expr-chart (or (:chart %) :cartesian)) datasets)
        groups (->> datasets
                    (map (fn [ds] {:m (:m ds) :alpha (:alpha ds) :ds ds}))
                    (group-by (fn [{:keys [m alpha]}] [(double m) (double alpha)]))
                    vals
                    (filter #(>= (count %) 2)))]
    (and (seq groups)
         (some (fn [grp]
                 (let [ms (keep identity (map #(conserved-mean-for-dataset ind %) (map :ds grp)))
                       n  (count ms)]
                   (when (>= n 2)
                     (let [gm   (/ (reduce + ms) n)
                           gstd (Math/sqrt
                                 (/ (reduce (fn [s v] (+ s (* (- v gm) (- v gm)))) 0.0 ms) n))]
                       (> (/ gstd (max 1e-10 (Math/abs gm))) 0.05)))))
               groups))))

(defn- expr-uses-position? [expr]
  (or (some #(expr-uses-symbol? expr %) '[qx qy])
      (some #(expr-uses-symbol? expr %) derived-vars)))

(defn- conserved-genome-valid?
  "A conserved-quantity genome must use at least one dynamic state variable
   (qx, qy, px, py) so trivially-constant or parameter-only expressions are rejected.
   We check the SIMPLIFIED form — GP sometimes hides px in canceling terms like
   (+ alpha px) - (c + px) which simplify to a constant.
   Also requires at least one position (qx/qy or derived r/r2/r3) AND one momentum
   (px/py) so qx-only flat functions cannot pass both CoV and ∇C·f checks.
   Requires IC dependence on the fixed reference scenarios (not just m or α)."
  [ind]
  (let [expr  (get ind conserved-expr-key)
        sexpr (when (some? expr) (simplify-expr expr))]
    (and (some? expr)
         (not (number? sexpr))
         (some #(expr-uses-symbol? sexpr %) '[qx qy px py r theta pr ptheta])
         (expr-uses-position? sexpr)
         (some #(expr-uses-symbol? sexpr %) '[px py pr ptheta])
         (try
           (ic-group-variation-ok? ind (reference-datasets))
           (catch Exception _ false)))))

(defn- single-law-valid? [ind]
  (case (law-kind ind)
    :analytical   (analytical-genome-valid? ind)
    :differential (symbols-covered-across-exprs? differential-expr-keys ind
                                                   required-differential-symbols)
    :conserved    (conserved-genome-valid? ind)
    false))

(defn genome-valid?
  "Analytical: each expr uses t; ICs and (m, α) appear across trajectory laws.
   Differential: state coords and (m, α) appear across rate laws.
   Conserved: expression uses at least one dynamic state variable.
   Composite: every law in :laws must be valid."
  [ind]
  (let [laws (individual-laws ind)]
    (and (seq laws)
         (every? #(single-law-valid? (law->legacy %)) laws))))

(defn- law-genome-key [law]
  (let [leg (law->legacy law)]
    (case (law-kind leg)
      :analytical   (into [:analytical]   (mapv leg analytical-expr-keys))
      :differential (into [:differential] (mapv leg differential-expr-keys))
      :conserved    [:conserved (get leg conserved-expr-key)]
      [(:strategy leg)])))

(defn individual-genome-key [ind]
  (into [:laws] (mapv law-genome-key (individual-laws ind))))

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

(defn- fitness-subsample-indices [n max-points]
  (cond
    (zero? n) []
    (<= n max-points) (vec (range n))
    (= max-points 1) [0]
    :else (vec (for [i (range max-points)]
                 (int (/ (* i (dec n)) (dec max-points)))))))

(defn- subsample-dataset
  [{:keys [data] :as dataset} max-points]
  (if (and max-points (pos? max-points) (> (count data) max-points))
    (assoc dataset :data (mapv #(nth data %) (fitness-subsample-indices (count data) max-points)))
    dataset))

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

(defn- conserved-probes-from-dataset [dataset]
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

(defn- conserved-probes-lite-from-dataset
  "One probe per scenario — enough to distinguish near-constant hacks, much cheaper than full sampling."
  [dataset]
  (when-let [s (first (:data dataset))]
    [{:qx (double (:qx s))
      :qy (double (:qy s))
      :px (double (:px s))
      :py (double (:py s))
      :m (double (:m dataset))
      :alpha (double (:alpha dataset))}]))

(defn build-behavior-probes
  "Probe tuples sampled from integrated scenarios — used to hash individuals by outputs, not syntax."
  [datasets]
  {:differential (vec (distinct (mapcat differential-probes-from-dataset datasets)))
   :analytical   (vec (distinct (mapcat analytical-probes-from-dataset datasets)))
   :conserved    (vec (distinct (mapcat conserved-probes-from-dataset datasets)))})

(defn build-behavior-probes-lite
  "Fewer probe points for per-generation escape/dedup (not final reporting)."
  [datasets]
  {:differential (vec (mapcat #(take 1 (differential-probes-from-dataset %)) datasets))
   :analytical   (vec (mapcat #(take 1 (analytical-probes-from-dataset %)) datasets))
   :conserved    (vec (mapcat conserved-probes-lite-from-dataset datasets))})

(defn- ^java.util.Map behavior-key-cache []
  (java.util.HashMap.))

(defn individual-behavior-key
  "Stable key from rounded model outputs on [[build-behavior-probes]]; nil if invalid or eval fails — use genome key then.
   Does not call [[genome-valid?]] — compile/eval failures return nil (expensive IC check skipped)."
  [ind {:keys [differential analytical conserved]}]
  (try
    (if (composite-individual? ind)
      (into [:composite]
            (keep #(individual-behavior-key (law->legacy %) {:differential differential
                                                             :analytical analytical
                                                             :conserved conserved})
                  (individual-laws ind)))
      (let [leg (law->legacy (first (individual-laws ind)))]
        (case (law-kind leg)
          :conserved
          (when (seq conserved)
            (let [c-fn (compile-conserved-fn (:c-expr leg))]
              [:conserved
               (vec (for [{:keys [qx qy px py m alpha]} conserved]
                      (quant-double (safe-double (c-fn qx qy px py m alpha)))))]))
          :differential
          (when (seq differential)
            (let [dqx (compile-rate-fn (:dqx-expr leg))
                  dqy (compile-rate-fn (:dqy-expr leg))
                  dpx (compile-rate-fn (:dpx-expr leg))
                  dpy (compile-rate-fn (:dpy-expr leg))]
              [:differential
               (vec
                (for [{:keys [qx qy px py m alpha]} differential]
                  [(quant-double (double (dqx qx qy px py m alpha)))
                   (quant-double (double (dqy qx qy px py m alpha)))
                   (quant-double (double (dpx qx qy px py m alpha)))
                   (quant-double (double (dpy qx qy px py m alpha)))]))]))
          :analytical
          (when (seq analytical)
            (let [qx-fn (compile-state-fn (:qx-expr leg))
                  qy-fn (compile-state-fn (:qy-expr leg))
                  px-fn (compile-state-fn (:px-expr leg))
                  py-fn (compile-state-fn (:py-expr leg))]
              [:analytical
               (vec
                (for [{:keys [t q0x q0y p0x p0y m alpha]} analytical]
                  [(quant-double (safe-double (qx-fn t q0x q0y p0x p0y m alpha)))
                   (quant-double (safe-double (qy-fn t q0x q0y p0x p0y m alpha)))
                   (quant-double (safe-double (px-fn t q0x q0y p0x p0y m alpha)))
                   (quant-double (safe-double (py-fn t q0x q0y p0x p0y m alpha)))]))]))
          nil)))
    (catch Exception _ nil)))

(defn- behavior-key-for
  "Behavior or genome key; cache is a mutable Map reused within one generation."
  [ind probes ^java.util.Map cache]
  (let [gk (individual-genome-key ind)]
    (if-some [hit (.get cache gk)]
      hit
      (let [k (or (individual-behavior-key ind probes) gk)]
        (.put cache gk k)
        k))))

(defn take-distinct-by-behavior
  "Keep top n ranked individuals with unique behavior key (fallback: genome key)."
  [n ranked probes & [cache]]
  (let [cache (or cache (behavior-key-cache))]
    (loop [seen #{} out [] xs (seq ranked)]
      (if (or (= (count out) n) (nil? xs))
        out
        (let [ind (first xs)
              k   (behavior-key-for ind probes cache)]
          (if (contains? seen k)
            (recur seen out (rest xs))
            (recur (conj seen k) (conj out ind) (rest xs))))))))

(defn- analytical-fitness-chart [ind _dataset]
  (infer-law-chart ind))

(defn- evaluate-analytical-predictions [ind dataset]
  (if-not (domain-applicable? ind dataset)
    {:strategy :analytical :n-horizons 0 :skipped? true
     :mse-q nil :mse-p nil :mse nil}
    (let [leg (if (per-regime-arm-fitness? ind)
                (analytical-law-regime-slice ind (dataset-regime-for-fitness dataset))
                ind)
          chart (analytical-fitness-chart leg dataset)]
      (if-not (coords/analytical-exprs-present? leg chart)
        {:strategy :analytical :n-horizons 0
         :mse-q Double/POSITIVE_INFINITY :mse-p Double/POSITIVE_INFINITY
         :mse Double/POSITIVE_INFINITY}
        (let [{:keys [data]} dataset
              times (all-horizon-times data)
              fns (compile-analytical-fns leg chart)
              ds  (assoc dataset :chart chart)
              {:keys [sq-q sq-p n]} (horizon-errors-analytical fns ds times)]
          (if (zero? n)
            {:strategy :analytical :n-horizons 0
             :mse-q Double/POSITIVE_INFINITY :mse-p Double/POSITIVE_INFINITY
             :mse Double/POSITIVE_INFINITY}
            (let [n (double n)]
              {:strategy :analytical :n-horizons (long n)
               :mse-q (/ sq-q n) :mse-p (/ sq-p n) :mse (/ (+ sq-q sq-p) n)})))))))

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

(defn- evaluate-conserved-predictions
  "Metrics for a conserved-quantity individual.
   :mse = max(squared CoV, squared relative RMS of |∇C·f|/|C|) — fails if either
   constancy or the conservation law is violated."
  [ind dataset]
  (try
    (let [chart (coords/effective-chart ind (:chart dataset))]
      (if-not (and (conserved-genome-valid? ind) (conserved-chart-compatible? ind chart))
        {:strategy :conserved :mse Double/POSITIVE_INFINITY :n-horizons 0}
        (let [ds     (assoc dataset :chart chart)
              c-fn   (compile-conserved-fn-for-chart chart (:c-expr ind))
              c-vals (conserved-vals-for-dataset c-fn ds)
              d-vals (conserved-dot-vals-for-dataset c-fn ds)
              fvals  (filterv #(Double/isFinite (safe-double %)) c-vals)
              n      (count fvals)]
          (if (< n 3)
            {:strategy :conserved :mse Double/POSITIVE_INFINITY :n-horizons 0}
            (let [mean     (/ (reduce + fvals) n)
                  variance (/ (reduce (fn [s v] (let [d (- v mean)] (+ s (* d d)))) 0.0 fvals) n)
                  cov-sq   (if (< (Math/abs mean) 1e-8) Double/POSITIVE_INFINITY
                             (/ variance (* mean mean)))
                  pairs    (filterv (fn [[c d]]
                                      (and (Double/isFinite (safe-double c))
                                           (Double/isFinite (safe-double d))))
                                  (map vector c-vals d-vals))
                  dot-sq   (if (empty? pairs) Double/POSITIVE_INFINITY
                             (/ (reduce (fn [s [c d]]
                                          (let [scale (max (Math/abs c) 1e-8)]
                                            (+ s (* (/ d scale) (/ d scale)))))
                                        0.0 pairs)
                                (count pairs)))]
              {:strategy :conserved
               :n-horizons n
               :mse (max cov-sq dot-sq)})))))
    (catch Exception _ {:strategy :conserved :mse Double/POSITIVE_INFINITY :n-horizons 0})))

(defn evaluate-predictions
  "Horizon error for an individual genome vs integrated trajectory (all t > 0).
   dataset is a scenario map with :data, :m, :alpha."
  [individual dataset]
  (case (:strategy individual)
    :analytical   (evaluate-analytical-predictions   individual dataset)
    :differential (evaluate-differential-predictions individual dataset)
    :conserved    (evaluate-conserved-predictions    individual dataset)
    {:mse Double/POSITIVE_INFINITY}))

(defn normalize-expr
  "GP trees must be plain lists for eval/print; mutate used to leave LazySeq subtrees."
  [expr]
  (cond
    (instance? clojure.lang.LazySeq expr) (normalize-expr (doall expr))
    (and (sequential? expr) (not (map? expr)))
    (apply list (map normalize-expr expr))
    :else expr))

(def ^:private latex-symbol-names
  '{t "t" m "m" alpha "\\alpha" energy "E"
    q0x "q_{0x}" q0y "q_{0y}" p0x "p_{0x}" p0y "p_{0y}"
    qx "q_x" qy "q_y" px "p_x" py "p_y"})

(def ^:private analytical-equation-specs
  [[:qx-expr :qx "q_x(t)" "q_x(t)"]
   [:qy-expr :qy "q_y(t)" "q_y(t)"]
   [:px-expr :px "p_x(t)" "p_x(t)"]
   [:py-expr :py "p_y(t)" "p_y(t)"]])

(def ^:private conserved-equation-specs
  [[:c-expr :c "C(q,p)" "C(\\mathbf{q},\\mathbf{p})"]])

(def ^:private differential-equation-specs
  [[:dqx-expr :dqx "dq_x/dt" "\\dot{q}_x"]
   [:dqy-expr :dqy "dq_y/dt" "\\dot{q}_y"]
   [:dpx-expr :dpx "dp_x/dt" "\\dot{p}_x"]
   [:dpy-expr :dpy "dp_y/dt" "\\dot{p}_y"]])

(defn- resolve-scenario
  "scenario — nil (first default), keyword :id, or full scenario map."
  [scenario]
  (cond
    (nil? scenario) (first default-scenarios)
    (keyword? scenario)
    (or (first (filter #(= (:id %) scenario) default-scenarios))
        (throw (ex-info "Unknown scenario id" {:scenario-id scenario
                                               :known (mapv :id default-scenarios)})))
    (and (map? scenario) (contains? scenario :q0x)) scenario
    :else (throw (ex-info "Expected scenario keyword or map with IC keys"
                          {:scenario scenario}))))

(defn- expr->math
  [expr & {:keys [latex?]}]
  (let [child #(expr->math % :latex? latex?)]
    (cond
      (number? expr) (str expr)
      (symbol? expr)
      (if latex?
        (get latex-symbol-names expr (name expr))
        (name expr))
      (sequential? expr)
      (let [[op & args] (normalize-expr expr)]
        (case op
          + (let [a (child (first args)) b (child (second args))]
              (if latex? (str a " + " b) (str "(" a " + " b ")")))
          - (let [a (child (first args)) b (child (second args))]
              (if latex? (str a " - " b) (str "(" a " - " b ")")))
          * (let [a (child (first args)) b (child (second args))]
              (if latex? (str a " \\cdot " b) (str "(" a " * " b ")")))
          e/square (let [a (child (first args))]
                     (if latex? (str a "^{2}") (str "square(" a ")")))
          e/sqrt (let [a (child (first args))]
                   (if latex? (str "\\sqrt{" a "}") (str "sqrt(" a ")")))
          e/sin (let [a (child (first args))]
                  (if latex? (str "\\sin\\left(" a "\\right)") (str "sin(" a ")")))
          e/cos (let [a (child (first args))]
                  (if latex? (str "\\cos\\left(" a "\\right)") (str "cos(" a ")")))
          e/sinh (let [a (child (first args))]
                   (if latex? (str "\\sinh\\left(" a "\\right)") (str "sinh(" a ")")))
          e/cosh (let [a (child (first args))]
                   (if latex? (str "\\cosh\\left(" a "\\right)") (str "cosh(" a ")")))
          e/mean-anomaly
          (let [a (child (first args)) b (child (second args)) c (child (nth args 2))]
            (if latex? (str "M(" a "," b "," c ")") (str "mean-anomaly(" a "," b "," c ")")))
          e/anomaly
          (let [a (child (first args)) b (child (second args)) c (child (nth args 2))]
            (if latex? (str "u(" a "," b "," c ")") (str "anomaly(" a "," b "," c ")")))
          e/periapsis-xp
          (let [s (clojure.string/join ", " (map child args))]
            (if latex? (str "x'(" s ")") (str "periapsis-xp(" s ")")))
          e/periapsis-yp
          (let [s (clojure.string/join ", " (map child args))]
            (if latex? (str "y'(" s ")") (str "periapsis-yp(" s ")")))
          e/lab-qx
          (let [s (clojure.string/join ", " (map child args))]
            (if latex? (str "q_x^{\\mathrm{lab}}(" s ")") (str "lab-qx(" s ")")))
          e/lab-qy
          (let [s (clojure.string/join ", " (map child args))]
            (if latex? (str "q_y^{\\mathrm{lab}}(" s ")") (str "lab-qy(" s ")")))
          e/div (let [a (child (first args)) b (child (second args))]
                  (if latex? (str "\\frac{" a "}{" b "}") (str "(" a " / " b ")")))
          e/if (let [[test then else] args
                     t (child test) a (child then) b (child else)]
                 (if latex?
                   (str "\\begin{cases} " a " & E<0 \\\\ " b " & E\\geq 0 \\end{cases}")
                   (str "if(" t "){" a "}{" b "}")))
          neg? (let [a (child (first args))]
                 (if latex? (str a " < 0") (str "neg?(" a ")")))
          (if latex?
            (str "\\mathrm{" (name op) "}\\left("
                 (clojure.string/join ", " (map child args)) "\\right)")
            (str "(" (name op) " " (clojure.string/join " " (map child args)) ")"))))
      :else (str expr))))

(defn- build-equation-map
  [individual specs format]
  (into {}
        (for [[expr-key var-key plain-lhs latex-lhs] specs
              :let [expr (normalize-expr (get individual expr-key))
                    lhs (if (= format :latex) latex-lhs plain-lhs)]]
          [var-key (str lhs " = " (expr->math expr :latex? (= format :latex)))])))

(defn- sample-analytical-at-t [individual dataset t]
  (let [{:keys [data m alpha]} dataset
        {:keys [q0x q0y p0x p0y]} (initial-ics data)
        md (double m)
        ad (double alpha)
        td (double t)
        fns {:qx (compile-state-fn (:qx-expr individual))
             :qy (compile-state-fn (:qy-expr individual))
             :px (compile-state-fn (:px-expr individual))
             :py (compile-state-fn (:py-expr individual))}
        predicted (into {}
                        (map (fn [[k f]]
                               [k (safe-double (f td q0x q0y p0x p0y md ad))])
                             fns))
        actual (when-let [s (state-at-t data t)]
                 (select-keys s [:qx :qy :px :py]))]
    {:t td
     :ics {:q0x q0x :q0y q0y :p0x p0x :p0y p0y :m md :alpha ad}
     :predicted predicted
     :actual actual}))

(defn- sample-differential-at-index [individual dataset idx]
  (let [{:keys [data m alpha]} dataset
        n (count data)
        i (long (max 0 (min (dec n) (or idx (quot n 2)))))
        s0 (nth data i)
        s1 (when (< i (dec n)) (nth data (inc i)))
        md (double m)
        ad (double alpha)
        fns {:dqx (compile-rate-fn (:dqx-expr individual))
             :dqy (compile-rate-fn (:dqy-expr individual))
             :dpx (compile-rate-fn (:dpx-expr individual))
             :dpy (compile-rate-fn (:dpy-expr individual))}
        rates (into {}
                    (map (fn [[k f]]
                           [k (safe-double (f (:qx s0) (:qy s0) (:px s0) (:py s0) md ad))])
                         fns))
        dt (when s1 (- (double (:t s1)) (double (:t s0))))
        predicted-next (when dt (symplectic-step-rates fns s0 dt m alpha))]
    {:index i
     :t (:t s0)
     :state (select-keys s0 [:qx :qy :px :py])
     :rates rates
     :dt dt
     :predicted-next (select-keys predicted-next [:qx :qy :px :py])
     :actual-next (when s1 (select-keys s1 [:qx :qy :px :py]))}))

(defn individual->equations
  "Readable equations for a genome plus optional numeric spot-check.

  individual — map with :laws or legacy :strategy and expression keys.
  scenario — nil, scenario :id keyword, or full scenario map (see resolve-scenario).
  opts:
    :format — :plain (default) or :latex
    :sample-t — analytical only: evaluate all four laws at this t
    :sample-index — differential only: trajectory index for state/rates (default: midpoint)

  Returns {:strategy :format :scenario-id :scenario-params :equations :metrics :sample}
  or {:laws [...]} for composite individuals."
  [individual & {:keys [scenario format sample-t sample-index]
                   :or {format :plain}}]
  (if (composite-individual? individual)
    {:laws (mapv #(individual->equations (law->legacy %) :scenario scenario :format format
                                         :sample-t sample-t :sample-index sample-index)
                  (individual-laws individual))}
    (let [ind (law->legacy (first (individual-laws individual)))
          sc (resolve-scenario scenario)
          dataset (scenario-data sc)
          specs (case (:strategy ind)
                  :analytical   analytical-equation-specs
                  :differential differential-equation-specs
                  :conserved    conserved-equation-specs
                  nil)
          _ (when (nil? specs)
              (throw (ex-info "Unknown individual strategy" {:strategy (:strategy ind)})))
          equations (build-equation-map ind specs format)
          metrics (evaluate-predictions ind dataset)
          sample (case (:strategy ind)
                   :analytical
                   (when sample-t
                     (sample-analytical-at-t ind dataset sample-t))
                   :differential
                   (sample-differential-at-index ind dataset sample-index)
                   :conserved nil)]
      {:strategy (:strategy ind)
       :format format
       :scenario-id (:id sc)
       :scenario-params (select-keys sc [:m :alpha :q0x :q0y :p0x :p0y :dt :steps])
       :equations equations
       :metrics metrics
       :sample sample})))

(def ops '[+ - * e/square e/sin e/cos e/sinh e/cosh e/div e/sqrt e/if])
;; Derived variables pre-computed by compile-state-fn from initial conditions.
;; r0/r02/r03/r05/r06 — initial radius powers (r₀¹,²,³,⁵,⁶)
;; omega             — circular angular velocity sqrt(α/m/r₀³); exact for circular orbit
;; omega-L           — actual initial angular velocity (q×p)/(m·r₀²); valid for any orbit
;; energy            — H(q₀,p₀); use with (e/if (neg? energy) bound-branch unbound-branch)
;; ecc, cos-om, sin-om, semi-a, semi-b, bh, n-mean — Kepler elements at t
;; kepler-u, kepler-F — eccentric/hyperbolic anomaly; kepler-xp/yp, kepler-vxp/vyp — periapsis frame
;; kepler-M0, kepler-M — mean anomaly at IC and at t (graded primitive pipeline)
(def analytical-derived-vars
  '[r0 r02 r03 r05 r06 omega omega-L energy
    ecc cos-om sin-om semi-a semi-b bh n-mean kepler-M0 kepler-M
    kepler-u kepler-F kepler-xp kepler-yp kepler-vxp kepler-vyp])
(def analytical-vars (vec (concat '[t] ic-vars param-vars analytical-derived-vars)))
(def differential-vars (vec (concat state-vars param-vars derived-vars)))

;; Variables for conserved-quantity expressions: current state, params, and
;; derived geometric variables (r, r2, r3 are pre-computed in compile-conserved-fn).
(def ^:private conserved-dynamic-vars
  (vec (concat state-vars '[r theta pr ptheta])))

(def ^:private conserved-momentum-vars '[px py pr ptheta])

(def conserved-vars (vec (concat conserved-dynamic-vars param-vars derived-vars)))

(def constants '[-1.0 -0.5 0.5 1.0 2.0 3.0])


(defn- random-atom [vars]
  (rand-nth (concat vars constants)))

(def unary-ops '#{e/square e/sin e/cos e/sinh e/cosh e/sqrt})
(def ternary-ops '#{e/if})

(defn- random-algebraic-expression
  [depth vars]
  (if (or (zero? depth) (< (rand) 0.3))
    (random-atom vars)
    (let [allowed (if (>= depth 2) ops (vec (remove #{'e/if} ops)))
          op (rand-nth allowed)]
      (cond
        (= op 'e/if)
        (list 'e/if '(neg? energy)
              (random-algebraic-expression (dec depth) vars)
              (random-algebraic-expression (dec depth) vars))

        (unary-ops op)
        (list op (random-algebraic-expression (dec depth) vars))

        :else
        (list op
              (random-algebraic-expression (dec depth) vars)
              (random-algebraic-expression (dec depth) vars))))))

(defn random-expression
  ([depth] (random-expression depth analytical-vars))
  ([depth vars]
   (if (and (pos? *primitive-tier*) (< (rand) primitive-use-rate))
     (let [p (prims/random-primitive-expr *primitive-tier* vars (min 2 depth))]
       (if (and p (prims/expr-primitive-valid? p *primitive-tier*))
         p
         (random-algebraic-expression depth vars)))
     (random-algebraic-expression depth vars))))

(defn- random-regime-branch-expr [depth vars]
  (list 'e/if '(neg? energy)
        (random-expression depth vars)
        (random-expression depth vars)))

(defn- random-valid-individual [genome-fn expr-keys required-syms valid?-fn]
  (loop [n 0]
    (let [ind (-> (genome-fn)
                  (ensure-symbol-coverage expr-keys required-syms))]
      (cond
        (valid?-fn ind) ind
        (> n 500) (if *both-regimes?*
                    (wrap-analytical-energy-branches
                     (into {:strategy :analytical :domain :any}
                           unbound-linear-analytical-exprs))
                    ind)
        :else (recur (inc n))))))

(defn random-analytical-individual
  ([] (random-analytical-individual (or *preferred-law-chart* :cartesian)))
  ([chart]
   (if (and (analytical-blocks-active?) (= chart :cartesian) (< (rand) 0.45))
     (catalog-law-as-individual)
     (let [chart    (keyword chart)
         keys     (coords/analytical-expr-keys-for-chart chart)
         vars     (vec (concat '[t]
                               (coords/ic-vars-for-chart chart)
                               param-vars
                               (if (= :polar chart)
                                 '[r02 r03 omega omega-L energy]
                                 analytical-derived-vars)))
         required (coords/required-analytical-symbols-for-chart chart)
         ensure-t (fn [ind]
                    (into ind
                          (map (fn [k]
                                 (let [expr (get ind k)]
                                   [k (if (expr-uses-t? expr) expr (list '+ 't expr))]))
                               keys)))]
     (random-valid-individual
      #(-> (into {:strategy :analytical
                  :domain (if *both-regimes?* :any
                              (if (< (rand) 0.35) (rand-nth [:bound :unbound]) :any))}
                 (map (fn [k]
                        [k (cond
                             (and *both-regimes?* (template-unbound-arms-active?))
                             (slot-with-bound-arm k (random-expression 4 vars))

                             (or *both-regimes?* (< (rand) 0.22))
                             (random-regime-branch-expr 3 vars)

                             :else
                             (random-expression 4 vars))])
                      keys))
           (ensure-symbol-coverage keys required)
           (cond-> (= chart :cartesian) ensure-analytical-uses-t)
           (cond-> (= chart :polar) ensure-t))
      keys
      required
      analytical-genome-valid?)))))

(defn random-differential-individual []
  (random-valid-individual
   #(into {:strategy :differential}
          (map (fn [k] [k (random-expression 4 differential-vars)])
               differential-expr-keys))
   differential-expr-keys
   required-differential-symbols
   #(genome-valid? %)))

;; ── Conserved Quantity Strategy ───────────────────────────────────────────────
;; Evolves a single expression C(qx, qy, px, py, m, α) that is constant along
;; any trajectory — i.e., a conservation law.  Fitness = 1/(CoV+1) where
;; CoV = std(C)/|mean(C)| along the trajectory.  Perfect conserved quantities
;; (energy, angular momentum) achieve fitness ≈ 1.0.
;;
;; conserved-expr-key, conserved-genome-valid?, and compile-conserved-fn are
;; defined earlier in the file so genome-valid? and evaluate-predictions
;; can reference them without forward declarations.

(defn- conserved-syntax-valid?
  "Cheap structural check — no per-scenario IC integration (used for escape immigrants)."
  [ind]
  (let [expr  (get ind conserved-expr-key)
        sexpr (when (some? expr) (simplify-expr expr))]
    (and (some? expr)
         (not (number? sexpr))
         (some #(expr-uses-symbol? sexpr %) '[qx qy px py r theta pr ptheta])
         (expr-uses-position? sexpr)
         (some #(expr-uses-symbol? sexpr %) '[px py pr ptheta]))))

(defn random-conserved-individual
  ([] (random-conserved-individual {}))
  ([{:keys [strict? max-tries] :or {strict? true max-tries 200}}]
   (let [valid? (if strict? conserved-genome-valid? conserved-syntax-valid?)
         max-tries (long max-tries)]
     (loop [n 0]
       (let [ind {:strategy :conserved
                  conserved-expr-key (random-expression 4 conserved-vars)}]
         (if (or (valid? ind) (>= n max-tries))
           ind
           (recur (inc n))))))))

(def ^:private conserved-rel-scale 12000.0)

(defn- fitness-from-relative-error
  "Sharper reward for near-zero relative error. Spurious near-constants with CoV or
   |∇C·f|/|C| ~0.01 score ~0.45; ~0.005 ~0.80; true invariants (≲1e-4) still ~1."
  [rel-metric]
  (/ 1.0 (+ 1.0 (* rel-metric rel-metric conserved-rel-scale))))

(defn- conserved-complexity-factor
  "Disfavor trig/sqrt trees — true Kepler invariants are compact algebra."
  [expr]
  (let [tree   (tree-seq sequential? seq (normalize-expr expr))
        trig-n (count (filter #(and (sequential? %) (#{'e/sin 'e/cos} (first %))) tree))
        sqrt-n (count (filter #(and (sequential? %) (= 'e/sqrt (first %))) tree))
        penalty (+ (* trig-n 0.25) (* sqrt-n 0.15))]
    (if (pos? penalty)
      (/ 1.0 (+ 1.0 penalty))
      1.0)))

(defn- fitness-from-conserved
  "fitness from CoV = std(C) / |mean(C)| along one trajectory."
  [vals]
  (let [finite (filterv #(Double/isFinite %) vals)
        n      (count finite)]
    (if (< n 3)
      0.0
      (let [mean     (/ (reduce + finite) n)
            variance (/ (reduce (fn [s v] (let [d (- v mean)] (+ s (* d d)))) 0.0 finite) n)
            std      (Math/sqrt variance)
            abs-mean (Math/abs mean)]
        (if (< abs-mean 1e-8)
          0.0
          (fitness-from-relative-error (/ std abs-mean)))))))

(defn- conserved-vals-for-dataset
  "Evaluate compiled conserved fn over every state in dataset; return vector of doubles."
  [c-fn {:keys [data m alpha chart] :as _ds}]
  (let [md (double m) ad (double alpha)]
    (if (= :polar chart)
      (mapv (fn [{:keys [r theta pr ptheta]}]
              (safe-double (c-fn (double r) (double theta) (double pr) (double ptheta) md ad)))
            data)
      (mapv (fn [{:keys [qx qy px py]}]
              (safe-double (c-fn (double qx) (double qy) (double px) (double py) md ad)))
            data))))

(def ^:private conserved-grad-eps 1e-6)

(defn- conserved-dot-at-state-polar
  [c-fn r theta pr ptheta m alpha]
  (let [md (double m) ad (double alpha)
        r (double r) theta (double theta) pr (double pr) ptheta (double ptheta)
        eps conserved-grad-eps
        c-at (fn [rv thv prv ptv]
               (safe-double (c-fn rv thv prv ptv md ad)))
        dc-dr (/ (- (c-at (+ r eps) theta pr ptheta) (c-at (- r eps) theta pr ptheta)) (* 2.0 eps))
        dc-dth (/ (- (c-at r (+ theta eps) pr ptheta) (c-at r (- theta eps) pr ptheta)) (* 2.0 eps))
        dc-dpr (/ (- (c-at r theta (+ pr eps) ptheta) (c-at r theta (- pr eps) ptheta)) (* 2.0 eps))
        dc-dpt (/ (- (c-at r theta pr (+ ptheta eps)) (c-at r theta pr (- ptheta eps))) (* 2.0 eps))
        {:keys [dr dtheta dpr dptheta]} (coords/polar-hamiltonian-deriv md ad {:r r :pr pr :ptheta ptheta})]
    (+ (* dc-dr dr) (* dc-dth dtheta) (* dc-dpr dpr) (* dc-dpt dptheta))))

(defn- conserved-dot-at-state-cartesian
  [c-fn qx qy px py m alpha]
  (let [md (double m) ad (double alpha)
        qx (double qx) qy (double qy) px (double px) py (double py)
        eps conserved-grad-eps
        c-at (fn [xq yq xp yp]
               (safe-double (c-fn xq yq xp yp md ad)))
        dc-dqx (/ (- (c-at (+ qx eps) qy px py) (c-at (- qx eps) qy px py)) (* 2.0 eps))
        dc-dqy (/ (- (c-at qx (+ qy eps) px py) (c-at qx (- qy eps) px py)) (* 2.0 eps))
        dc-dpx (/ (- (c-at qx qy (+ px eps) py) (c-at qx qy (- px eps) py)) (* 2.0 eps))
        dc-dpy (/ (- (c-at qx qy px (+ py eps)) (c-at qx qy px (- py eps))) (* 2.0 eps))
        {:keys [dqx dqy dpx dpy]} (grav2d-deriv md ad {:qx qx :qy qy :px px :py py})]
    (+ (* dc-dqx dqx) (* dc-dqy dqy) (* dc-dpx dpx) (* dc-dpy dpy))))

(defn- conserved-dot-at-state
  ([c-fn qx qy px py m alpha]
   (conserved-dot-at-state-cartesian c-fn qx qy px py m alpha))
  ([c-fn state m alpha chart]
   (if (= :polar chart)
     (let [{:keys [r theta pr ptheta]} state]
       (conserved-dot-at-state-polar c-fn r theta pr ptheta m alpha))
     (let [{:keys [qx qy px py]} state]
       (conserved-dot-at-state-cartesian c-fn qx qy px py m alpha)))))

(defn- conserved-dot-vals-for-dataset
  "Per-state |dC/dt| = |∇C·f| along a trajectory."
  [c-fn {:keys [data m alpha chart] :as _ds}]
  (let [md (double m) ad (double alpha)
        chart (or chart :cartesian)]
    (mapv (fn [state]
            (let [dot (conserved-dot-at-state c-fn state md ad chart)]
              (if (Double/isFinite dot) (Math/abs dot) Double/NaN)))
          data)))

(defn- fitness-from-conservation-law
  "fitness from RMS(|∇C·f| / |C|) — true invariants have ∇C·f = 0."
  [c-vals dot-vals]
  (let [pairs  (filterv (fn [[c d]] (and (Double/isFinite c) (Double/isFinite d)))
                        (map vector c-vals dot-vals))
        n      (count pairs)]
    (if (< n 3)
      0.0
      (let [sq-rel (reduce (fn [s [c d]]
                             (let [scale (max (Math/abs c) 1e-8)]
                               (+ s (* (/ d scale) (/ d scale)))))
                           0.0 pairs)
            rms-rel (Math/sqrt (/ sq-rel n))]
        (fitness-from-relative-error rms-rel)))))

(defn- conserved-trajectory-fitness
  "Combined per-trajectory fitness: must be both constant (low CoV) AND a true
   invariant (∇C·f ≈ 0). Returns the minimum of the two component scores."
  [c-fn dataset]
  (let [c-vals (conserved-vals-for-dataset c-fn dataset)
        d-vals (conserved-dot-vals-for-dataset c-fn dataset)]
    (min (fitness-from-conserved c-vals)
         (fitness-from-conservation-law c-vals d-vals))))

(defn- random-phase-state
  "One point in (qx, qy, px, py, m, α) inside [[default-scenario-bounds]]."
  [& {:keys [bounds rng]}]
  (let [bounds (or bounds default-scenario-bounds)
        {:keys [m alpha q0x q0y p0x p0y r-min]} bounds
        r-min (or r-min 1.0)]
    (loop [attempt 0]
      (let [qx (uniform-sample rng (first q0x) (second q0x))
            qy (uniform-sample rng (first q0y) (second q0y))
            r  (Math/sqrt (+ (* qx qx) (* qy qy)))]
        (if (>= r r-min)
          {:qx (double qx) :qy (double qy)
           :px (uniform-sample rng (first p0x) (second p0x))
           :py (uniform-sample rng (first p0y) (second p0y))
           :m  (uniform-sample rng (first m) (second m))
           :alpha (uniform-sample rng (first alpha) (second alpha))}
          (if (> attempt 200)
            (throw (ex-info "Could not sample valid phase state |q| >= r-min"
                            {:bounds bounds :r-min r-min}))
            (recur (inc attempt))))))))

(defn phase-states-for-fitness-context
  "Materialize phase-space sample points for one evaluation batch (:de-driven)."
  [ctx & {:keys [generation]}]
  (let [n (:phase-samples ctx 48)
        seed (when (:seed ctx) (+ (long (:seed ctx)) (long (or generation 0))))
        rng  (when seed (java.util.Random. seed))]
    (mapv (fn [_] (random-phase-state :rng rng)) (range n))))

(def ^:private analytical-ode-check-times
  "Short times after each IC for analytical ODE-residual probes."
  [0.0 0.005 0.01 0.02 0.04])

(defn- analytical-state-at
  [fns t q0x q0y p0x p0y m alpha]
  (let [td (double t) md (double m) ad (double alpha)
        {qx-fn :qx qy-fn :qy px-fn :px py-fn :py} fns]
    {:qx (safe-double (qx-fn td q0x q0y p0x p0y md ad))
     :qy (safe-double (qy-fn td q0x q0y p0x p0y md ad))
     :px (safe-double (px-fn td q0x q0y p0x p0y md ad))
     :py (safe-double (py-fn td q0x q0y p0x p0y md ad))}))

(defn- analytical-ode-residual-sq
  "Squared ODE residual ‖(dq/dt, dp/dt)_num − f(q,p)‖² at one (IC, t) probe."
  [fns {:keys [t q0x q0y p0x p0y m alpha]}]
  (let [dt   1e-4
        td   (double t)
        md   (double m)
        ad   (double alpha)
        s    (analytical-state-at fns td q0x q0y p0x p0y md ad)
        {:keys [qx qy px py]} s
        deriv (grav2d-deriv md ad s)
        [dqx-num dqy-num dpx-num dpy-num]
        (if (< td (* 2.0 dt))
          (let [s1 (analytical-state-at fns (+ td dt) q0x q0y p0x p0y md ad)]
            [(/ (- (:qx s1) qx) dt)
             (/ (- (:qy s1) qy) dt)
             (/ (- (:px s1) px) dt)
             (/ (- (:py s1) py) dt)])
          (let [sm (analytical-state-at fns (- td dt) q0x q0y p0x p0y md ad)
                sp (analytical-state-at fns (+ td dt) q0x q0y p0x p0y md ad)]
            [(/ (- (:qx sp) (:qx sm)) (* 2.0 dt))
             (/ (- (:qy sp) (:qy sm)) (* 2.0 dt))
             (/ (- (:px sp) (:px sm)) (* 2.0 dt))
             (/ (- (:py sp) (:py sm)) (* 2.0 dt))]))
        D (max 1.0 (Math/sqrt (+ (* qx qx) (* qy qy))))
        scale-p (max 1.0 (Math/sqrt (+ (* px px) (* py py))))]
    (+ (e/square (/ (- dqx-num (:dqx deriv)) D))
       (e/square (/ (- dqy-num (:dqy deriv)) D))
       (e/square (/ (- dpx-num (:dpx deriv)) scale-p))
       (e/square (/ (- dpy-num (:dpy deriv)) scale-p)))))

(defn- fitness-from-mean-sq [mean-sq]
  (if (or (not (Double/isFinite mean-sq)) (neg? mean-sq))
    0.0
    (/ 1.0 (+ 1.0 (Math/sqrt mean-sq)))))

(defn- calculate-analytical-equation-fitness
  "Score analytical trajectories by ODE residual d/dt q,p ≈ f(q,p) on random ICs."
  [ind phase-states]
  (try
    (when (analytical-genome-valid? ind)
      (let [fns {:qx (compile-state-fn (:qx-expr ind))
                 :qy (compile-state-fn (:qy-expr ind))
                 :px (compile-state-fn (:px-expr ind))
                 :py (compile-state-fn (:py-expr ind))}
            probes (vec (mapcat
                         (fn [{:keys [qx qy px py m alpha]}]
                           (map (fn [t]
                                  {:t t :q0x qx :q0y qy :p0x px :p0y py :m m :alpha alpha})
                                analytical-ode-check-times))
                         phase-states))
            sqs    (filterv #(and (Double/isFinite %) (not (Double/isNaN %)))
                            (map #(analytical-ode-residual-sq fns %) probes))]
        (when (seq sqs)
          (fitness-from-mean-sq (/ (reduce + sqs) (count sqs))))))
    (catch Exception _ 0)))

(defn- calculate-conserved-equation-fitness
  "Score conserved C by Poisson condition ∇C·f ≈ 0 at phase-space samples (true ODE in f)."
  [ind phase-states]
  (try
    (when (conserved-genome-valid? ind)
      (let [c-fn   (compile-conserved-fn (:c-expr ind))
            c-vals (mapv (fn [{:keys [qx qy px py m alpha]}]
                           (safe-double (c-fn (double qx) (double qy)
                                              (double px) (double py)
                                              (double m) (double alpha))))
                         phase-states)
            finite (filterv #(Double/isFinite %) c-vals)
            n      (count finite)]
        (when (>= n 3)
          (let [mean   (/ (reduce + finite) n)
                variance (/ (reduce (fn [s v] (let [d (- v mean)] (+ s (* d d)))) 0.0 finite) n)
                cstd   (Math/sqrt variance)
                cov-ok? (> (/ cstd (max 1e-10 (Math/abs mean))) 0.05)
                dot-vals (mapv (fn [{:keys [qx qy px py m alpha]}]
                                 (conserved-dot-at-state c-fn qx qy px py m alpha))
                               phase-states)
                law-fit (fitness-from-conservation-law c-vals dot-vals)]
            (when (and cov-ok? (pos? law-fit))
              (* (conserved-complexity-factor (:c-expr ind))
                 law-fit))))))
    (catch Exception _ 0)))

(defn calculate-conserved-fitness
  "Score a conserved-quantity individual on one trajectory dataset.
   Requires low CoV along the orbit AND ∇C·f ≈ 0 (true invariant)."
  [ind dataset]
  (try
    (if (and (conserved-genome-valid? ind)
             (conserved-chart-compatible? ind (coords/effective-chart ind (:chart dataset))))
      (let [chart (coords/effective-chart ind (:chart dataset))
            c-fn  (compile-conserved-fn-for-chart chart (:c-expr ind))]
        (conserved-trajectory-fitness c-fn (assoc dataset :chart chart)))
      0.0)
    (catch Exception _ 0.0)))

(def ^:dynamic *fast-conserved-fitness?* false)
(def ^:dynamic *fast-differential-fitness?* false)
(def ^:dynamic *fitness-max-states* nil)
(def ^:dynamic *fitness-timeout-ms* nil)

(defn- calculate-conserved-trajectory-de-fitness
  "DE-driven conserved: C must be constant in time along integrated orbits
   (low CoV AND ∇C·f ≈ 0 on each reference trajectory)."
  [ind datasets]
  (try
    (when (and (seq datasets) (conserved-genome-valid? ind))
      (let [datasets (if-let [ms *fitness-max-states*]
                       (mapv #(subsample-dataset % ms) datasets)
                       datasets)
            fits   (mapv (fn [ds]
                           (let [chart (coords/effective-chart ind (:chart ds))]
                             (when (conserved-chart-compatible? ind chart)
                               (let [c-fn (compile-conserved-fn-for-chart chart (:c-expr ind))]
                                 (conserved-trajectory-fitness c-fn (assoc ds :chart chart))))))
                         datasets)
            fits   (filterv some? fits)]
        (when (and (seq fits) (every? pos? fits))
          (* (conserved-complexity-factor (:c-expr ind))
             (apply min fits)))))
    (catch Exception _ nil)))

(defn- calculate-conserved-de-driven-fitness
  "Conserved laws are independent invariants, not solutions of the motion DE.
   Score only temporal constancy and ∇C·f ≈ 0 along integrated physical trajectories."
  [ind _phase-states datasets]
  (or (calculate-conserved-trajectory-de-fitness ind datasets) 0.0))

(defn- call-with-timeout [ms label f]
  (if (and ms (pos? ms))
    (let [fut (future (try (f) (catch Exception _ 0.0)))
          res (deref fut ms ::timeout)]
      (when (= res ::timeout)
        (future-cancel fut true)
        (println (str "    fitness TIMEOUT after " ms "ms"
                      (when label (str " (" label ")"))))
        (flush))
      (if (= res ::timeout) 0.0 res))
    (f)))

(defn- calculate-conserved-fitness-scenarios
  "Full conserved-quantity fitness:
   1. Compile expression once.
   2. Per-trajectory: low CoV AND |∇C·f| ≈ 0 (Poisson / Noether condition).
   3. IC-group check on fixed reference scenarios (reject parameter-only hacks).
   4. Global cross-scenario variation on the training/eval batch.
   Returns worst per-trajectory combined fitness after all checks pass."
  [ind datasets]
  (try
    (when (if *fast-conserved-fitness?*
            (conserved-syntax-valid? ind)
            (conserved-genome-valid? ind))
      (let [datasets (if-let [ms *fitness-max-states*]
                       (mapv #(subsample-dataset % ms) datasets)
                       datasets)
            per-traj (mapv (fn [ds]
                             (let [chart (coords/effective-chart ind (:chart ds))]
                               (when (conserved-chart-compatible? ind chart)
                                 (let [ds'  (assoc ds :chart chart)
                                       c-fn (compile-conserved-fn-for-chart chart (:c-expr ind))
                                       vals (conserved-vals-for-dataset c-fn ds')
                                       finite (filterv #(Double/isFinite %) vals)
                                       n      (count finite)
                                       mean   (when (>= n 1) (/ (reduce + finite) n))]
                                   {:vals finite :n n :mean mean :m (:m ds) :alpha (:alpha ds)
                                    :ds ds' :chart chart :c-fn c-fn}))))
                           datasets)
            per-traj (filterv some? per-traj)
            traj-means (keep :mean per-traj)
            ic-ok?   (if *fast-conserved-fitness?*
                       true
                       (ic-group-variation-ok? ind (reference-datasets)))
            global-ok? (let [fmeans (filterv some? traj-means)
                             n      (count fmeans)]
                         (and (>= n 2)
                              (let [cm   (/ (reduce + fmeans) n)
                                    cstd (Math/sqrt
                                          (/ (reduce (fn [s v] (+ s (* (- v cm) (- v cm)))) 0.0 fmeans) n))]
                                (> (/ cstd (max 1e-10 (Math/abs cm))) 0.05))))]
        (if-not (and ic-ok? global-ok?)
          0.0
          (* (conserved-complexity-factor (:c-expr ind))
             (apply min (map (fn [{:keys [c-fn ds]}]
                               (conserved-trajectory-fitness c-fn ds))
                             per-traj))))))
    (catch Exception _ 0)))

(def ^:dynamic *strategy-filter*
  "When set to a keyword (:analytical, :differential, :conserved), random-individual
   generates only that strategy.  nil (default) = uniform mix of all three."
  nil)

(def ^:dynamic *fast-immigrants?*
  "When true, conserved immigrants skip expensive IC integration (fitness filters invalid)."
  false)

(def ^:dynamic *fast-breeding?*
  "When true, skip block decomposition / behavior dedup during breeding (escape mode)."
  false)

(def ^:dynamic *de-driven-search?*
  "When true, immigrants and default strategy mix exclude :differential (redundant with known ODE)."
  false)

(def default-mcts-mutate-rate 0.2)
(def default-mcts-mutate-simulations 48)

(def ^:dynamic *mcts-mutate?*
  "When true, analytical mutate-individual may run MCTS slot repair instead of GP tree walk."
  false)

(def ^:dynamic *mcts-mutate-rate*
  "Probability of attempting MCTS repair on each analytical mutation (not restart/crossover)."
  default-mcts-mutate-rate)

(def ^:dynamic *mcts-mutate-simulations*
  "MCTS simulations per mutation attempt."
  default-mcts-mutate-simulations)

(def ^:dynamic *mcts-mutate-datasets*
  "Scenario datasets for MCTS mutation fitness (set during evolve-generation)."
  nil)

(defn- random-law [kind]
  (law-from-legacy
   (case kind
     :analytical   (random-analytical-individual)
     :differential (random-differential-individual)
     :conserved    (if *fast-immigrants?*
                      (random-conserved-individual {:strict? false :max-tries 25})
                      (random-conserved-individual)))))

(defn- random-composite-de-individual []
  (laws-individual [(random-analytical-individual)
                    (random-law :conserved)]))

(defn random-individual []
  (cond
    ;; DE-driven: one genome with analytical trajectories + a conserved invariant.
    (and *de-driven-search?*
         (or (nil? *strategy-filter*)
             (= *strategy-filter* :conserved)))
    (random-composite-de-individual)

    (and *de-driven-search?* (= *strategy-filter* :analytical))
    (laws-individual [(random-law :analytical)])

    *strategy-filter*
    (laws-individual [(random-law *strategy-filter*)])

    :else
    (laws-individual
     [(random-law
       (let [r (rand)]
         (cond
           (< r 0.34) :analytical
           (< r 0.67) :differential
           :else :conserved)))])))

(defn- escape-deep-reseed-population
  "Full population for [escape+]: cycle validated catalog laws (strict template), else random."
  [n]
  (if (and (= *strategy-filter* :analytical) (analytical-blocks-active?))
    (let [entries (cond-> (vec blocks/validated-catalog)
                    (blocks/kepler-conic-valid?) (conj (blocks/kepler-conic-entry)))
          inds (mapv #(laws-individual [(law-from-legacy (strict-catalog-law-from-entry %))])
                     entries)]
      (if (seq inds)
        (vec (take n (cycle inds)))
        (vec (repeatedly n random-individual))))
    (vec (repeatedly n random-individual))))

;; ── Physics seeds ─────────────────────────────────────────────────────────────
;; Hand-crafted individuals encoding known correct or near-correct physics.
;; Injected into the initial population when --seed is passed so the GP can
;; refine from a good starting point rather than discovering 1/r³ from scratch.

(def physics-seeds
  "Seed individuals that express known gravitational physics (or close variants).
   All expressions use (* -1.0 alpha) for negation since the GP only has binary
   ops; unary (- alpha) would be mangled by simplify-expr."
  [;; Exact Newtonian gravity: dq/dt = p/m, dp/dt = -α·q/r³
   {:strategy :differential
    :dqx-expr '(e/div px m)
    :dqy-expr '(e/div py m)
    :dpx-expr '(e/div (* (* -1.0 alpha) qx) r3)
    :dpy-expr '(e/div (* (* -1.0 alpha) qy) r3)}
   ;; Same but r3 written as r·r2 — different subtree structure for crossover diversity
   {:strategy :differential
    :dqx-expr '(e/div px m)
    :dqy-expr '(e/div py m)
    :dpx-expr '(e/div (* (* -1.0 alpha) qx) (* r r2))
    :dpy-expr '(e/div (* (* -1.0 alpha) qy) (* r r2))}
   ;; Perturbed coefficient -0.9 — lets evolution tune from slightly-off value
   {:strategy :differential
    :dqx-expr '(e/div px m)
    :dqy-expr '(e/div py m)
    :dpx-expr '(e/div (* (* -0.9 alpha) qx) r3)
    :dpy-expr '(e/div (* (* -0.9 alpha) qy) r3)}
   ;; Perturbed coefficient -1.1
   {:strategy :differential
    :dqx-expr '(e/div px m)
    :dqy-expr '(e/div py m)
    :dpx-expr '(e/div (* (* -1.1 alpha) qx) r3)
    :dpy-expr '(e/div (* (* -1.1 alpha) qy) r3)}
   ;; Alternative form: force = -alpha*(q/r)/r2 — different factorisation
   {:strategy :differential
    :dqx-expr '(e/div px m)
    :dqy-expr '(e/div py m)
    :dpx-expr '(e/div (* (* -1.0 alpha) (e/div qx r)) r2)
    :dpy-expr '(e/div (* (* -1.0 alpha) (e/div qy r)) r2)}

   ;; ── Analytical (short-horizon Taylor) seeds ───────────────────────────
   ;; First-order: q(t) ≈ q₀ + (p₀/m)·t,  p(t) ≈ p₀ + F(q₀)·t
   ;; r0/r02/r03 are pre-computed by compile-state-fn from initial position.
   {:strategy :analytical
    :domain :bound
    :qx-expr '(+ q0x (* (e/div p0x m) t))
    :qy-expr '(+ q0y (* (e/div p0y m) t))
    :px-expr '(+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
    :py-expr '(+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))}

   ;; Second-order: adds ½·F·t²/m correction to positions
   ;; q(t) ≈ q₀ + (p₀/m)t − (α·q₀)/(2m·r₀³)·t²
   {:strategy :analytical
    :domain :bound
    :qx-expr '(+ q0x (+ (* (e/div p0x m) t) (* (e/div (* (* -0.5 alpha) q0x) (* m r03)) (* t t))))
    :qy-expr '(+ q0y (+ (* (e/div p0y m) t) (* (e/div (* (* -0.5 alpha) q0y) (* m r03)) (* t t))))
    :px-expr '(+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
    :py-expr '(+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))}

   ;; Circular rotation seed (exact for circular orbit, approximation for ellipse/hyperbola):
   ;; q(t) = R(ωt)·q₀  where ω = sqrt(α/(m·r₀³)), pre-computed as `omega`.
   ;; p(t) = m·ω·R(ωt+π/2)·q₀ = m·ω·(−sin,cos)·q₀
   {:strategy :analytical
    :domain :bound
    :qx-expr '(- (* q0x (e/cos (* omega t))) (* q0y (e/sin (* omega t))))
    :qy-expr '(+ (* q0x (e/sin (* omega t))) (* q0y (e/cos (* omega t))))
    :px-expr '(* (* -1.0 (* m omega)) (+ (* q0x (e/sin (* omega t))) (* q0y (e/cos (* omega t)))))
    :py-expr '(* (* m omega) (- (* q0x (e/cos (* omega t))) (* q0y (e/sin (* omega t)))))}

   ;; Angular-momentum rotation seed: uses the ACTUAL initial angular velocity Ω_L = L/(m·r₀²)
   ;; where L = q₀×p₀ = q0x·p0y − q0y·p0x.  Better than circular-omega for eccentric orbits.
   ;; Momentum follows from d/dt of the rotated position.
   {:strategy :analytical
    :domain :bound
    :qx-expr '(- (* q0x (e/cos (* omega-L t))) (* q0y (e/sin (* omega-L t))))
    :qy-expr '(+ (* q0x (e/sin (* omega-L t))) (* q0y (e/cos (* omega-L t))))
    :px-expr '(* (* -1.0 (* m omega-L)) (+ (* q0x (e/sin (* omega-L t))) (* q0y (e/cos (* omega-L t)))))
    :py-expr '(* (* m omega-L) (- (* q0x (e/cos (* omega-L t))) (* q0y (e/sin (* omega-L t)))))}

   ;; Taylor2 with 2nd-order momentum correction:
   ;; q(t) ≈ q₀ + (p₀/m)t − ½·(α·q₀/m)·t²/r₀³
   ;; p(t) ≈ p₀ + F₀·t + ½·(dF/dt)|₀·t²
   ;; where dFx/dt|₀ = α·[(-1/r₀³ + 3q₀x²/r₀⁵)·p₀x + 3q₀x·q₀y·p₀y/r₀⁵] / m
   {:strategy :analytical
    :domain :bound
    :qx-expr '(+ q0x (+ (* (e/div p0x m) t) (* (e/div (* (* -0.5 alpha) q0x) (* m r03)) (* t t))))
    :qy-expr '(+ q0y (+ (* (e/div p0y m) t) (* (e/div (* (* -0.5 alpha) q0y) (* m r03)) (* t t))))
    :px-expr '(+ p0x (+ (* (e/div (* (* -1.0 alpha) q0x) r03) t)
                        (* 0.5
                           (* (e/div alpha (* m r05))
                              (* (+ (* (- (* 3.0 (* q0x q0x)) r02) p0x)
                                    (* (* 3.0 q0x) (* q0y p0y)))
                                 t))
                           t)))
    :py-expr '(+ p0y (+ (* (e/div (* (* -1.0 alpha) q0y) r03) t)
                        (* 0.5
                           (* (e/div alpha (* m r05))
                              (* (+ (* (- (* 3.0 (* q0y q0y)) r02) p0y)
                                    (* (* 3.0 q0y) (* q0x p0x)))
                                 t))
                           t)))}

   ;; Third-order momentum + second-order position seed:
   ;; q(t) as before (2nd-order Taylor); p(t) adds 3rd-order correction.
   ;; d³px/dt³|₀ ≈ α²·q₀x / (m·r₀⁶)   (circular-orbit approximation for the curvature term)
   ;; This fixes the remaining 16% error in px for the circle scenario.
   {:strategy :analytical
    :domain :bound
    :qx-expr '(+ q0x (+ (* (e/div p0x m) t) (* (e/div (* (* -0.5 alpha) q0x) (* m r03)) (* t t))))
    :qy-expr '(+ q0y (+ (* (e/div p0y m) t) (* (e/div (* (* -0.5 alpha) q0y) (* m r03)) (* t t))))
    :px-expr '(+ p0x (+ (* (e/div (* (* -1.0 alpha) q0x) r03) t)
                        (+ (* 0.5 (* (e/div alpha (* m r05))
                                     (* (+ (* (- (* 3.0 (* q0x q0x)) r02) p0x)
                                           (* (* 3.0 q0x) (* q0y p0y)))
                                        t))
                               t)
                           (* (e/div (* alpha (* alpha q0x)) (* 6.0 (* m r06)))
                              (* t (* t t))))))
    :py-expr '(+ p0y (+ (* (e/div (* (* -1.0 alpha) q0y) r03) t)
                        (+ (* 0.5 (* (e/div alpha (* m r05))
                                     (* (+ (* (- (* 3.0 (* q0y q0y)) r02) p0y)
                                           (* (* 3.0 q0y) (* q0x p0x)))
                                        t))
                               t)
                           (* (e/div (* alpha (* alpha q0y)) (* 6.0 (* m r06)))
                              (* t (* t t))))))}

   ;; Regime branching: Taylor on bound orbits, same linearization on unbound (no :domain tag).
   {:strategy :analytical
    :domain :any
    :qx-expr '(e/if (neg? energy)
                 (+ q0x (* (e/div p0x m) t))
                 (+ q0x (* (e/div p0x m) t)))
    :qy-expr '(e/if (neg? energy)
                 (+ q0y (* (e/div p0y m) t))
                 (+ q0y (* (e/div p0y m) t)))
    :px-expr '(e/if (neg? energy)
                 (+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
                 (+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t)))
    :py-expr '(e/if (neg? energy)
                 (+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))
                 (+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t)))}

   ;; ── Conserved quantity seeds ──────────────────────────────────────────
   ;; Angular momentum: L = qx·py − qy·px  (always conserved in central force)
   {:strategy :conserved
    :c-expr '(- (* qx py) (* qy px))}
   ;; Energy: H = (px²+py²)/(2m) − α/r  (Hamiltonian, conserved by definition)
   {:strategy :conserved
    :c-expr '(- (e/div (+ (* px px) (* py py)) (* 2.0 m)) (e/div alpha r))}
   ;; Kinetic energy alone — wrong but nearby in search space
   {:strategy :conserved
    :c-expr '(e/div (+ (* px px) (* py py)) (* 2.0 m))}])

(defn- both-regimes-analytical-seeds []
  (mapv #(-> %
               (cond-> (not (analytical-branches-on-energy? %))
                 (wrap-analytical-energy-branches))
               template-unbound-arms-law)
        (into (vec (filterv #(= (:strategy %) :analytical) physics-seeds))
              (when *analytical-blocks?*
                (blocks/catalog-laws)))))

(def de-driven-composite-seeds
  "Bundled analytical q(t),p(t) + Hamiltonian energy — one individual, two laws."
  [{:laws [{:kind :analytical
            :domain :bound
            :r-expr '(+ r0 (* (e/div pr0 m) t))
            :theta-expr '(+ theta0 (* (e/div ptheta0 (* m r02)) t))
            :pr-expr '(+ pr0 (* (e/div (* (* -1.0 alpha) ptheta0 ptheta0) (* m r03 r02)) t))
            :ptheta-expr 'ptheta0}
           {:kind :conserved
            :c-expr '(- (e/div (+ (* pr pr) (* ptheta ptheta (* r r))) (* 2.0 m)) (e/div alpha r))}]}])


(defn- char-scale
  "Characteristic length scale of a dataset: RMS of |q| over the trajectory.
   Used to make RMSE dimensionless so fitness = 1 / (RMSE/D + 1) ∈ (0, 1]."
  [data]
  (let [n (count data)]
    (if (zero? n)
      1.0
      (Math/sqrt (/ (reduce + (map (fn [{:keys [qx qy]}]
                                     (+ (* (double qx) (double qx))
                                        (* (double qy) (double qy))))
                                   data))
                    (double n))))))

(defn- fitness-from-errors
  "fitness = 1 / (RMSE/D + 1)  where D is the characteristic position scale.
   Result is dimensionless and bounded in (0, 1]; returns 0 for degenerate inputs."
  [{:keys [sq-q sq-p n]} D]
  (if (or (zero? n) (<= (double D) 0.0))
    0.0
    (let [n (double n)
          sum-sq (+ sq-q sq-p)]
      (if (or (Double/isNaN sum-sq) (Double/isInfinite sum-sq))
        0.0
        (let [rmse (Math/sqrt (/ sum-sq n))]
          (if (Double/isNaN rmse)
            0.0
            (/ 1.0 (+ (/ rmse (double D)) 1.0))))))))

(declare calculate-analytical-slice-fitness-at-horizon)

(defn calculate-analytical-fitness-at-horizon
  "Single-scenario analytical fitness at one horizon fraction, or dual training horizons
   (20% + full orbit) when `horizon-frac` is omitted."
  ([ind dataset]
   (apply min
          (mapv #(calculate-analytical-fitness-at-horizon ind dataset %)
                analytical-training-horizon-fracs)))
  ([ind dataset horizon-frac]
   (try
     (cond
       (not (analytical-genome-valid? ind))
       0.0

       (per-regime-arm-fitness? ind)
       (calculate-analytical-slice-fitness-at-horizon
        (analytical-law-regime-slice ind (dataset-regime-for-fitness dataset))
        dataset horizon-frac)

       :else
       (let [chart (analytical-fitness-chart ind dataset)]
         (if (coords/analytical-exprs-present? ind chart)
           (let [{:keys [data]} dataset
                 D    (char-scale data)
                 fns  (compile-analytical-fns ind chart)
                 ds   (assoc dataset :chart chart)
                 errors (horizon-errors-analytical fns ds (horizon-fraction-times data horizon-frac))]
             (fitness-from-errors errors D))
           0.0)))
     (catch Exception _ 0.0))))

(defn calculate-analytical-fitness [ind dataset]
  (calculate-analytical-fitness-at-horizon ind dataset))

(defn adversarial-repair-target
  "For an analytical individual, return the weakest scenario and expression slot to patch.
   {:dataset :expr-key :scenario-fitness :mse-q :mse-p :scenario-id}"
  [ind datasets]
  (let [leg (first-law-legacy ind)]
    (when (= :analytical (:strategy leg))
      (let [applicable (filterv #(domain-applicable? leg %) datasets)]
        (when (seq applicable)
          (let [{:keys [dataset scenario-fitness]}
                (apply min-key :scenario-fitness
                       (mapv (fn [ds]
                               {:dataset ds
                                :scenario-fitness
                                (calculate-analytical-fitness-at-horizon leg ds)})
                             applicable))
                metrics (evaluate-predictions leg dataset)
                finite? (fn [x] (and (number? x) (Double/isFinite x) (pos? x)))
                mse-q (:mse-q metrics)
                mse-p (:mse-p metrics)
                expr-key (cond
                           (and (finite? mse-q) (finite? mse-p))
                           (if (>= mse-q mse-p)
                             (rand-nth [:qx-expr :qy-expr])
                             (rand-nth [:px-expr :py-expr]))
                           (finite? mse-q) (rand-nth [:qx-expr :qy-expr])
                           (finite? mse-p) (rand-nth [:px-expr :py-expr])
                           :else (rand-nth analytical-expr-keys))]
            {:dataset dataset
             :expr-key expr-key
             :scenario-fitness scenario-fitness
             :mse-q mse-q
             :mse-p mse-p
             :scenario-id (:id (:scenario dataset))}))))))

(defn apply-analytical-repair
  "Insert a repaired analytical law into an individual (re-wraps energy branches when needed)."
  [ind repaired-leg]
  (let [leg (if (and *both-regimes?* (not (analytical-strict-energy-branches? repaired-leg)))
              (wrap-analytical-energy-branches repaired-leg)
              repaired-leg)
        law (law-from-legacy (template-unbound-arms-law leg))]
    (if (seq (:laws ind))
      (assoc-in ind [:laws 0] law)
      (law->legacy law))))

(defn calculate-differential-fitness [ind dataset]
  (try
    (when (genome-valid? ind)
      (let [{:keys [data]} dataset
            D   (char-scale data)
            fns {:dqx (compile-rate-fn (:dqx-expr ind))
                 :dqy (compile-rate-fn (:dqy-expr ind))
                 :dpx (compile-rate-fn (:dpx-expr ind))
                 :dpy (compile-rate-fn (:dpy-expr ind))}
            dataset (if-let [ms *fitness-max-states*]
                      (subsample-dataset dataset ms)
                      dataset)
            one-step (one-step-errors-differential fns dataset)
            errors   (if *fast-differential-fitness?*
                       one-step
                       (let [rollout (rollout-errors-differential fns dataset)]
                         {:sq-q (+ (:sq-q one-step) (:sq-q rollout))
                          :sq-p (+ (:sq-p one-step) (:sq-p rollout))
                          :n    (+ (:n one-step)    (:n rollout))}))]
        (fitness-from-errors errors D)))
    (catch Exception _ 0)))

(defn calculate-fitness
  "dataset is a scenario map with :data, :m, :alpha."
  [individual dataset]
  (or (case (primary-strategy-label individual)
        :analytical   (calculate-analytical-fitness   (law->legacy (first (individual-laws individual))) dataset)
        :differential (calculate-differential-fitness (law->legacy (first (individual-laws individual))) dataset)
        :conserved    (calculate-conserved-fitness    (law->legacy (first (individual-laws individual))) dataset)
        :composite
        (let [fits (mapv #(calculate-fitness (law->legacy %) dataset)
                         (individual-laws individual))]
          (if (seq fits) (apply min fits) 0.0))
        nil)
      0))

(defn aggregate-scenario-fitness
  "Combine per-scenario fitness values: :min (default) or :percentile."
  [fits & {:keys [aggregate percentile] :or {aggregate :min percentile 0.1}}]
  (if (empty? fits)
    0.0
    (let [sorted (sort fits)]
      (case (keyword aggregate)
        :min (double (first sorted))
        :percentile
        (let [n (count sorted)
              idx (min (dec n) (max 0 (int (Math/floor (* (double percentile) (dec n))))))]
          (double (nth sorted idx)))
        (throw (ex-info "Unknown fitness aggregate" {:aggregate aggregate}))))))

(defn- calculate-analytical-slice-fitness-at-horizon
  "Fitness for a regime-sliced law (arms only, no outer e/if) on one scenario."
  [ind dataset horizon-frac]
  (try
    (let [chart (analytical-fitness-chart ind dataset)]
      (if-not (analytical-slice-valid? ind chart)
        0.0
        (if-not (coords/analytical-exprs-present? ind chart)
          0.0
          (let [{:keys [data]} dataset
                D    (char-scale data)
                fns  (compile-analytical-fns ind chart)
                ds   (assoc dataset :chart chart)
                errors (horizon-errors-analytical fns ds (horizon-fraction-times data horizon-frac))]
            (fitness-from-errors errors D)))))
    (catch Exception _ 0.0)))

(defn- calculate-analytical-slice-fitness
  "Regime-sliced fitness: min over 20% and full-orbit horizons on one scenario."
  [ind dataset]
  (apply min
         (mapv #(calculate-analytical-slice-fitness-at-horizon ind dataset %)
               analytical-training-horizon-fracs)))

(defn- calculate-analytical-per-regime-fitness
  "Both-regimes: bound arms on bound scenarios, unbound arms on unbound; fitness = min(regimes).
   Each scenario scored at min(20% horizon, full orbit)."
  [ind datasets & {:keys [aggregate percentile] :or {aggregate :min percentile 0.1}}]
  (let [bound-ds   (filterv #(= :bound (dataset-regime-for-fitness %)) datasets)
        unbound-ds (filterv #(= :unbound (dataset-regime-for-fitness %)) datasets)
        bound-ind  (analytical-law-regime-slice ind :bound)
        unbound-ind (analytical-law-regime-slice ind :unbound)
        bound-fit  (when (seq bound-ds)
                     (aggregate-scenario-fitness
                      (mapv #(calculate-analytical-slice-fitness bound-ind %)
                            bound-ds)
                      :aggregate aggregate :percentile percentile))
        unbound-fit (when (seq unbound-ds)
                      (aggregate-scenario-fitness
                       (mapv #(calculate-analytical-slice-fitness unbound-ind %)
                             unbound-ds)
                       :aggregate aggregate :percentile percentile))
        regime-fits (vec (remove nil? [bound-fit unbound-fit]))]
    (if (seq regime-fits)
      (aggregate-scenario-fitness regime-fits :aggregate :min)
      0.0)))

(defn- calculate-analytical-de-driven-fitness
  "DE-driven analytical: 20% + full-orbit trajectory fit on integrated reference orbits.
   With both-regimes strict branches: bound arms on bound scenarios, unbound on unbound."
  [ind datasets & {:keys [aggregate percentile] :or {aggregate :min percentile 0.1}}]
  (try
    (when (and (seq datasets) (analytical-genome-valid? ind))
      (let [datasets (datasets-for-analytical-law ind datasets)]
        (if (empty? datasets)
          0.0
          (if (per-regime-arm-fitness? ind)
            (calculate-analytical-per-regime-fitness
             ind datasets :aggregate aggregate :percentile percentile)
            (aggregate-scenario-fitness
             (mapv (fn [ds] (calculate-analytical-fitness ind ds)) datasets)
             :aggregate aggregate :percentile percentile)))))
    (catch Exception _ nil)))

(defn- calculate-de-driven-fitness
  "Fitness from the environment physics — not naked trajectory matching in data-driven mode.
   Analytical: 20% + full-orbit fit on integrated orbits. Conserved: invariance along orbits only.
   Differential rate laws score 0 (redundant with the known equations of motion)."
  [ind _phase-states datasets & {:keys [aggregate percentile] :or {aggregate :min percentile 0.1}}]
  (case (:strategy ind)
    :analytical   (or (calculate-analytical-de-driven-fitness ind datasets
                                                               :aggregate aggregate
                                                               :percentile percentile) 0.0)
    :conserved    (or (calculate-conserved-de-driven-fitness ind nil datasets) 0.0)
    :differential 0.0
    0.0))

(defn- calculate-single-fitness-scenarios
  "Fitness for one law (legacy single-strategy map)."
  [individual datasets & {:keys [aggregate percentile evaluation phase-states]
                          :or {aggregate :min percentile 0.1 evaluation :data-driven}}]
  (if (= (keyword evaluation) :de-driven)
    (or (calculate-de-driven-fitness individual phase-states datasets
                                     :aggregate aggregate :percentile percentile) 0.0)
    (if (= (:strategy individual) :conserved)
      (or (calculate-conserved-fitness-scenarios individual datasets) 0.0)
      (let [datasets (if (= (:strategy individual) :analytical)
                       (datasets-for-analytical-law individual datasets)
                       datasets)]
        (if (empty? datasets)
          0.0
          (aggregate-scenario-fitness
           (mapv #(calculate-fitness individual %) datasets)
           :aggregate aggregate
           :percentile percentile))))))

(defn calculate-fitness-scenarios
  "Robust fitness over scenario datasets: :min (worst case) or :percentile (e.g. p10).

  When :evaluation is :de-driven, scores analytical laws by min(20% horizon, full-orbit)
  trajectory fit on integrated reference orbits; conserved laws by invariance along those
  orbits only (not as
  DE solutions). Differential always scores 0.

  Composite individuals (:laws with multiple entries) take the worst law fitness.

  For :conserved in :data-driven mode, uses a combined per-trajectory + cross-scenario check
  so trivially-constant expressions (e.g. cos(cos(qy/qy)) = const) score 0.

  Optional :aggregate and :percentile (see [[aggregate-scenario-fitness]])."
  [individual datasets & {:keys [aggregate percentile evaluation phase-states]
                          :or {aggregate :min percentile 0.1 evaluation :data-driven}}]
  (let [laws (individual-laws individual)]
    (if (<= (count laws) 1)
      (calculate-single-fitness-scenarios
       (if (seq laws) (law->legacy (first laws)) individual)
       datasets
       :aggregate aggregate :percentile percentile :evaluation evaluation
       :phase-states phase-states)
      (aggregate-scenario-fitness
       (mapv #(calculate-single-fitness-scenarios
               (law->legacy %)
               datasets
               :aggregate aggregate :percentile percentile :evaluation evaluation
               :phase-states phase-states)
             laws)
       :aggregate :min))))

(defn- expr-subtrees [expr]
  (let [expr (normalize-expr expr)]
    (if (sequential? expr)
      (cons expr (mapcat expr-subtrees (rest expr)))
      [expr])))

(defn- slotted-subtrees [ind]
  (mapcat
   (fn [law]
     (let [leg  (law->legacy law)
           keys (case (law-kind leg)
                  :analytical   analytical-expr-keys
                  :differential differential-expr-keys
                  :conserved    [conserved-expr-key]
                  [])]
       (for [slot keys, st (expr-subtrees (get leg slot))]
         {:slot slot :law (law-kind leg) :motif st})))
   (individual-laws ind)))

(def ^:private known-physics-motifs
  "Canonical 2D gravity subtrees → short gloss (exact structural match after normalize-expr)."
  {'(e/div px m) "Hamilton: q̇_x = p_x/m"
   '(e/div py m) "Hamilton: q̇_y = p_y/m"
   '(e/div p0x m) "IC form: q̇_x ≈ p_{0x}/m at t₀"
   '(e/div p0y m) "IC form: q̇_y ≈ p_{0y}/m at t₀"})

(defn- motif-tags [motif]
  (cond-> #{}
    (expr-uses-symbol? motif 'm) (conj :uses-m)
    (expr-uses-symbol? motif 'alpha) (conj :uses-alpha)
    (expr-uses-symbol? motif 't) (conj :uses-t)
    (some #(expr-uses-symbol? motif %) ic-vars) (conj :uses-ic)
    (some #(expr-uses-symbol? motif %) state-vars) (conj :uses-state)))

(defn- motif-gloss [motif slot]
  (or (get known-physics-motifs (normalize-expr motif))
      (cond
        (and (#{:dqx-expr :dqy-expr} slot) (seq (motif-tags motif)))
        (str "rate slot " (name slot) " — "
             (clojure.string/join ", " (map name (motif-tags motif))))
        (and (#{:qx-expr :qy-expr :px-expr :py-expr} slot) (expr-uses-symbol? motif 't))
        "trajectory uses time"
        :else nil)))

(defn population->motif-report
  "Mine recurring subtrees (micro-heuristics) from a ranked or unranked population.

  population — vector of individuals (e.g. from [[load-population]]).
  datasets — scenario maps with :data (defaults to [[scenarios->datasets]] on [[default-scenarios]]).
  opts:
    :elite-frac — top fraction by fitness to mine (default 0.5)
    :min-share   — minimum fraction of elite carrying motif (default 0.2)
    :top         — max motifs to return (default 15)

  Returns {:n-population :n-elite :elite-min-fitness :motifs [...] :physics-hits [...]}."
  [population & {:keys [datasets elite-frac min-share top]
                 :or {elite-frac 0.5 min-share 0.2 top 15}}]
  (let [datasets (or datasets (scenarios->datasets default-scenarios))
        ranked (->> population
                    (map (fn [ind]
                           (assoc ind :fitness (or (:fitness ind)
                                                   (calculate-fitness-scenarios ind datasets)))))
                    (sort-by :fitness >))
        n-pop (count ranked)
        n-elite (max 1 (int (Math/ceil (* elite-frac n-pop))))
        elite (take n-elite ranked)
        min-fit (if (seq elite) (:fitness (last elite)) 0.0)
        by-motif (reduce
                  (fn [acc ind]
                    (let [fit (:fitness ind)
                          ;; one hit per individual per [slot motif], not per subtree occurrence
                          unique (reduce (fn [m {:keys [slot motif]}]
                                           (let [sig (pr-str (normalize-expr motif))]
                                             (assoc m [slot sig] (normalize-expr motif))))
                                         {}
                                         (slotted-subtrees ind))]
                      (reduce
                       (fn [a [[slot sig] motif]]
                         (update a [slot sig]
                                   (fn [v]
                                     (let [base (or v {:slot slot
                                                       :signature sig
                                                       :motif motif
                                                       :count 0
                                                       :fitness-sum 0.0})]
                                       (-> base
                                           (update :count inc)
                                           (update :fitness-sum + fit))))))
                       acc
                       unique)))
                  {}
                  elite)
        enriched (for [[_ v] by-motif
                       :let [cnt (:count v)
                             share (/ (double cnt) n-elite)
                             avg (/ (:fitness-sum v) (double cnt))]
                       :when (>= share min-share)]
                   (let [motif (:motif v)]
                     (assoc v
                            :share share
                            :avg-fitness avg
                            :tags (motif-tags motif)
                            :gloss (motif-gloss motif (:slot v)))))
        sorted (take top (sort-by (juxt :share :avg-fitness) #(compare %2 %1) enriched))
        physics-hits (vec (keep #(when (get known-physics-motifs (:motif %))
                                  (assoc (select-keys % [:slot :motif :signature :share])
                                         :gloss (get known-physics-motifs (:motif %))))
                            enriched))]
    {:n-population n-pop
     :n-elite n-elite
     :elite-min-fitness min-fit
     :motifs sorted
     :physics-hits physics-hits}))

(defn- perturb-constants
  "Walk an expression tree and jitter every numeric constant by ±20%.
   Leaves zero as zero to avoid introducing small biases."
  [expr]
  (normalize-expr
   (cond
     (and (number? expr) (not (zero? expr)))
     (* expr (+ 1.0 (* 0.4 (- (rand) 0.5))))   ; ×(0.8 … 1.2)
     (coll? expr)
     (apply list (first expr) (map perturb-constants (rest expr)))
     :else expr)))

(defn mutate
  ([expr] (mutate expr analytical-vars))
  ([expr vars]
   (normalize-expr
    (let [r (rand)]
      (cond
        ;; 15% — replace whole subtree with a fresh random expression
        (< r 0.15) (random-expression 2 vars)
        ;; 10% — jitter all numeric constants in the expression
        (< r 0.25) (perturb-constants expr)
        ;; 75% — recurse into children, mutating each independently
        (coll? expr)
        (let [op (first expr)]
          (cons op (mapv #(mutate % vars) (rest expr))))
        :else expr)))))

;; ── Functional blocks (modular GP) ───────────────────────────────────────────
;; Each expression slot decomposes into composable blocks (+/−/* terms or one subtree).
;; Crossover swaps blocks between genomes — e.g. kinetic block + potential block → H.

(defn- individual-expr-keys [ind-or-law]
  (let [k (or (:kind ind-or-law) (:strategy ind-or-law))]
    (case k
      :analytical   (coords/analytical-expr-keys-for-chart (infer-law-chart ind-or-law))
      :differential differential-expr-keys
      :conserved    [conserved-expr-key]
      [])))

(defn- decompose-expr-to-blocks
  "Split an expression into a composable block representation."
  [expr]
  (let [e (normalize-expr expr)]
    (cond
      (not (sequential? e)) {:compose :mono :blocks [e]}
      (= '+ (first e)) {:compose :+ :blocks (vec (rest e))}
      (and (= '- (first e)) (= 1 (count (rest e))))
      {:compose :neg :blocks [(second e)]}
      (and (= '- (first e)) (>= (count (rest e)) 2))
      {:compose :- :blocks [(second e) (nth e 2)]}
      (= '* (first e)) {:compose :* :blocks (vec (rest e))}
      :else {:compose :mono :blocks [e]})))

(defn- compose-blocks-to-expr
  [{:keys [compose blocks]}]
  (let [blocks (vec blocks)]
    (when (seq blocks)
      (normalize-expr
       (case compose
         :mono (first blocks)
         :neg  (list '- (first blocks))
         :-   (list '- (first blocks) (second blocks))
         :+    (apply list '+ blocks)
         :*    (apply list '* blocks)
         (first blocks))))))

(def ^:private max-blocks-per-slot 6)

(defn- attach-block-structure [ind]
  (assoc ind :expr-blocks
         (into {}
               (for [k (individual-expr-keys ind)]
                 [k (decompose-expr-to-blocks (get ind k))]))))

(defn- materialize-exprs-from-blocks [ind]
  (reduce (fn [acc [slot blks]]
            (if-let [expr (compose-blocks-to-expr blks)]
              (assoc acc slot expr)
              acc))
          ind
          (:expr-blocks ind)))

(defn sync-block-genome
  "Keep :expr-blocks aligned with expression keys.
   When :from-blocks? true, materialize expressions from blocks first (after crossover)."
  [ind & {:keys [from-blocks?]}]
  (if (and from-blocks? (:expr-blocks ind))
    (-> ind materialize-exprs-from-blocks attach-block-structure)
    (attach-block-structure ind)))

(defn- sync-law [law]
  (law-from-legacy (sync-block-genome (law->legacy law))))

(defn sync-individual [ind]
  (if (seq (:laws ind))
    (assoc ind :laws (mapv sync-law (:laws ind)))
    (sync-block-genome ind)))

(defn- block-crossover-blocks
  "Graft block(s) from donor into recipient — no external catalog."
  [recv donor]
  (let [recv (update recv :blocks vec)
        rb   (vec (:blocks recv))
        db   (vec (:blocks donor))
        graft (when (seq db) (rand-nth db))]
    (if (nil? graft)
      recv
      (let [op (rand-nth [:replace :append :adopt-compose])]
        (case op
          :replace
          (if (seq rb)
            (assoc recv :blocks (assoc rb (rand-int (count rb)) graft))
            {:compose :mono :blocks [graft]})

          :append
          (cond
            (= :mono (:compose recv))
            {:compose :+ :blocks [(first rb) graft]}
            (and (= :+ (:compose recv)) (< (count rb) max-blocks-per-slot))
            (update recv :blocks conj graft)
            (seq rb)
            (assoc recv :blocks (assoc rb (rand-int (count rb)) graft))
            :else {:compose :mono :blocks [graft]})

          :adopt-compose
          (if (and (#{:+ :- :*} (:compose donor)) (seq db))
            (assoc donor :blocks (vec (take max-blocks-per-slot db)))
            recv))))))

(defn- mutate-block-repr [block-repr vars]
  (let [blocks (vec (:blocks block-repr))]
    (if (empty? blocks)
      block-repr
      (let [i (rand-int (count blocks))]
        (update block-repr :blocks assoc i (mutate (nth blocks i) vars))))))

(defn- mutate-expr-via-blocks [expr vars _slot]
  (-> expr
      decompose-expr-to-blocks
      (mutate-block-repr vars)
      compose-blocks-to-expr))

(defn individual-block-summary
  "Block decomposition per slot (dev / REPL)."
  [ind]
  (let [ind (sync-block-genome ind)]
    (into {} (for [k (individual-expr-keys ind)] [k (get-in ind [:expr-blocks k])]))))

;; ── Symbolic block abstraction + hypothesis mutations ─────────────────────────
;; Classify blocks by role (kinetic / potential / angular / junk) and apply
;; structural edits using blocks already present in the genome — recompose,
;; drop junk, swap in sibling blocks. No physics catalog or template injection.

(def ^:dynamic *guess-mutations?* true)
;; Set during long stagnation below target fitness — more random restarts, fewer clones.
(def ^:dynamic *stagnation-escape?* false)

(defn- expr-tree-size [expr]
  (count (tree-seq sequential? seq (normalize-expr expr))))

(defn- expr-uses-op? [expr op]
  (cond
    (= expr op) true
    (coll? expr) (boolean (some #(expr-uses-op? % op) expr))
    :else false))

(defn- angular-like? [expr]
  (let [s (normalize-expr (simplify-expr expr))]
    (or (= s '(- (* qx py) (* qy px)))
        (and (expr-uses-symbol? s 'qx) (expr-uses-symbol? s 'py)
             (expr-uses-symbol? s 'qy) (expr-uses-symbol? s 'px)
             (not (expr-uses-op? s 'e/sin))
             (not (expr-uses-op? s 'e/cos))))))

(defn- kinetic-like? [expr]
  (and (expr-uses-symbol? expr 'm)
       (or (expr-uses-symbol? expr 'px) (expr-uses-symbol? expr 'py))
       (or (expr-uses-op? expr 'e/square)
           (some #(and (sequential? %) (= 'e/square (first %)))
                 (tree-seq sequential? seq expr)))))

(defn- potential-like? [expr]
  (and (expr-uses-symbol? expr 'alpha)
       (expr-uses-position? expr)
       (not (expr-uses-symbol? expr 'px))
       (not (expr-uses-symbol? expr 'py))))

(defn- block-junk-score [expr]
  (+ (* 2.0 (count (filter #(and (sequential? %) (#{'e/sin 'e/cos 'e/sinh 'e/cosh} (first %)))
                          (tree-seq sequential? seq expr))))
     (* 1.0 (count (filter #(and (sequential? %) (= 'e/sqrt (first %)))
                          (tree-seq sequential? seq expr))))
     (if (and (expr-uses-symbol? expr 'm)
              (> (expr-tree-size expr) 8)
              (not (or (kinetic-like? expr) (potential-like? expr) (angular-like? expr))))
       2.0
       0.0)))

(defn abstract-block
  "Symbolic summary of one functional block (not a full seed — local classification)."
  [expr]
  (let [expr (normalize-expr expr)]
    {:expr expr
     :size (expr-tree-size expr)
     :tags (motif-tags expr)
     :kinetic? (kinetic-like? expr)
     :potential? (potential-like? expr)
     :angular? (angular-like? expr)
     :trig? (or (expr-uses-op? expr 'e/sin) (expr-uses-op? expr 'e/cos)
                (expr-uses-op? expr 'e/sinh) (expr-uses-op? expr 'e/cosh))
     :junk-score (block-junk-score expr)}))

(defn abstract-slot-blocks
  "Symbolic profile of all blocks in one expression slot."
  [block-repr]
  (let [blocks (mapv abstract-block (:blocks block-repr))]
    {:compose (:compose block-repr)
     :blocks blocks
     :has-kinetic (some :kinetic? blocks)
     :has-potential (some :potential? blocks)
     :has-angular (some :angular? blocks)
     :has-trig (some :trig? blocks)
     :max-junk (if (seq blocks) (apply max (map :junk-score blocks)) 0.0)
     :junkiest-idx (when (seq blocks)
                     (first (apply max-key (fn [[_i b]] (:junk-score b))
                                         (map-indexed vector blocks))))}))

(defn block-abstraction-summary
  "Per-slot symbolic block profiles for an individual (dev / REPL)."
  [ind]
  (let [ind (sync-block-genome ind)]
    (into {}
          (for [k (individual-expr-keys ind)]
            [k (abstract-slot-blocks (get-in ind [:expr-blocks k]))]))))

(defn- replace-block-at-index [block-repr idx new-expr]
  (update block-repr :blocks assoc idx (normalize-expr new-expr)))

(defn- drop-block-at-index [block-repr idx]
  (let [nb (vec (concat (subvec (vec (:blocks block-repr)) 0 idx)
                        (subvec (vec (:blocks block-repr)) (inc idx))))]
    (when (seq nb)
      (if (= 1 (count nb))
        {:compose :mono :blocks nb}
        (assoc block-repr :blocks nb)))))

(defn- conserved-guess-candidates [block-repr]
  (let [prof       (abstract-slot-blocks block-repr)
        blocks     (:blocks prof)
        k-expr     (some :expr (filter :kinetic? blocks))
        p-expr     (some :expr (filter :potential? blocks))
        good-exprs (mapv :expr (filter #(or (:kinetic? %) (:potential? %) (:angular? %)) blocks))
        cands      (atom [])]
    (when (and k-expr p-expr)
      (swap! cands conj {:compose :- :blocks [k-expr p-expr]}))
    (when (and (>= (count blocks) 2) (= :+ (:compose block-repr)))
      (swap! cands conj (assoc block-repr :compose :-)))
    (when (and (>= (count blocks) 2) (= :- (:compose block-repr)))
      (swap! cands conj (assoc block-repr :compose :+)))
    (when-let [ji (:junkiest-idx prof)]
      (when-let [dropped (drop-block-at-index block-repr ji)]
        (swap! cands conj dropped))
      (doseq [g (distinct good-exprs)]
        (swap! cands conj (replace-block-at-index block-repr ji g))))
    (let [non-trig (vec (map :expr (filter #(not (:trig? %)) blocks)))]
      (when (>= (count non-trig) 2)
        (swap! cands conj {:compose :- :blocks [(first non-trig) (second non-trig)]})
        (swap! cands conj {:compose :+ :blocks [(first non-trig) (second non-trig)]})))
    (distinct (filter some? @cands))))

(defn- differential-guess-candidates [block-repr]
  (let [prof       (abstract-slot-blocks block-repr)
        blocks     (:blocks prof)
        good-exprs (mapv :expr (remove #(> (:junk-score %) 1.0) blocks))
        cands      (atom [])]
    (when (and (>= (count blocks) 2) (= :+ (:compose block-repr)))
      (swap! cands conj (assoc block-repr :compose :-)))
    (when-let [ji (:junkiest-idx prof)]
      (when-let [dropped (drop-block-at-index block-repr ji)]
        (swap! cands conj dropped))
      (doseq [g (distinct good-exprs)]
        (swap! cands conj (replace-block-at-index block-repr ji g))))
    (distinct (filter some? @cands))))

(defn- slot-guess-candidates [slot block-repr]
  (let [catalog-cands
        (when (analytical-blocks-active?)
          (vec (for [entry (blocks/entries-for-slot slot)]
                 (decompose-expr-to-blocks (blocks/slot-expr entry slot)))))]
    (into (or catalog-cands [])
          (case slot
            :c-expr (conserved-guess-candidates block-repr)
            (:dqx-expr :dqy-expr :dpx-expr :dpy-expr) (differential-guess-candidates block-repr)
            []))))

(defn guess-mutate-individual
  "Apply one symbolically guessed block edit (hypothesis mutation), not a random tree walk."
  [ind]
  (let [ind (sync-block-genome ind)
        ks  (individual-expr-keys ind)]
    (if (empty? ks)
      ind
      (let [slot (rand-nth ks)
            blks (get-in ind [:expr-blocks slot])
            cands (seq (slot-guess-candidates slot blks))]
        (if cands
          (let [full (filter #(#{:- :mono} (:compose %)) (vec cands))
                pick (if (and (seq full) (< (rand) 0.55))
                       (rand-nth full)
                       (rand-nth cands))]
            (sync-block-genome
             (assoc-in ind [:expr-blocks slot] pick)
             :from-blocks? true))
          ind)))))

(defn- random-block-mutate-individual [ind]
  (let [ind (sync-block-genome ind)
        mutate-one (fn [acc k vars]
                     (if (< (rand) 0.7)
                       (update-in acc [:expr-blocks k] mutate-block-repr vars)
                       (assoc-in acc [:expr-blocks k]
                                 (decompose-expr-to-blocks
                                  (mutate (get acc k) vars)))))]
    (sync-block-genome
     (case (:strategy ind)
       :analytical
       (reduce #(mutate-one %1 %2 analytical-vars) ind analytical-expr-keys)
       :differential
       (reduce #(mutate-one %1 %2 differential-vars) ind differential-expr-keys)
       :conserved
       (mutate-one ind conserved-expr-key conserved-vars)
       ind)
     :from-blocks? true)))

(defn- validate-and-repair-law
  [ind]
  (let [strategy (:strategy ind)
        [ks vars] (case strategy
                    :analytical   [analytical-expr-keys   analytical-vars]
                    :differential [differential-expr-keys differential-vars]
                    :conserved    [[conserved-expr-key]   conserved-vars]
                    [nil nil])
        ind' (if ks
               (reduce (fn [acc k]
                         (let [raw (get acc k)
                               expr (cond
                                      (expr-too-large? raw) (random-expression 3 vars)
                                      :else (let [s (simplify-if-small raw)]
                                              (if (number? s)
                                                (random-expression 3 vars)
                                                s)))]
                           (assoc acc k expr)))
                       ind ks)
               ind)
        repaired (case strategy
                   :analytical   (-> ind'
                                     (ensure-symbol-coverage ks required-analytical-symbols)
                                     ensure-analytical-uses-t)
                   :differential (ensure-symbol-coverage ind' ks required-differential-symbols)
                   :conserved    (if (conserved-genome-valid? ind')
                                   ind'
                                   (random-conserved-individual))
                   ind')
        ind''  (if *fast-breeding?* repaired (sync-block-genome repaired))]
    (template-unbound-arms-law ind'')))

(defn- validate-and-repair
  "Simplify, detect degenerate (constant) sub-expressions, and ensure required
   symbols survive after mutation/crossover."
  [ind]
  (if (seq (:laws ind))
    (laws-individual (mapv #(law-from-legacy (validate-and-repair-law (law->legacy %)))
                            (:laws ind)))
    (validate-and-repair-law ind)))

(defn- mutate-bound-arm-at-slot [ind k]
  (let [bound (regime-arm-expr (normalize-expr (get ind k)) :bound)
        mutated (mutate bound analytical-vars)]
    (assoc ind k (slot-with-bound-arm k mutated))))

(defn- guess-mutate-bound-arm [ind k]
  (let [bound (regime-arm-expr (normalize-expr (get ind k)) :bound)
        temp  (-> ind
                  (assoc k bound)
                  sync-block-genome)
        blks  (get-in temp [:expr-blocks k])
        cands (seq (slot-guess-candidates k blks))]
    (if cands
      (let [pick (rand-nth cands)
            synced (sync-block-genome (assoc-in temp [:expr-blocks k] pick) :from-blocks? true)]
        (assoc ind k (slot-with-bound-arm k (get synced k))))
      ind)))

(defn- mutate-bound-slot [ind k]
  (if (and *guess-mutations?* (< (rand) 0.45))
    (let [g (guess-mutate-bound-arm ind k)]
      (if (not= g ind) g (mutate-bound-arm-at-slot ind k)))
    (mutate-bound-arm-at-slot ind k)))

(defn- pair-mutate-analytical-individual [ind]
  (let [pair (rand-nth bound-arm-pairs)]
    (-> (cond
          (and (analytical-blocks-active?) (< (rand) 0.3))
          (apply-catalog-pair-to-law ind)

          (and (template-conic-unbound-active?) (< (rand) 0.25))
          (apply-catalog-unbound-pair-to-law ind)

          :else
          (reduce mutate-bound-slot ind pair))
        template-unbound-arms-law
        validate-and-repair-law)))

(defn- mutate-analytical-domain [law]
  (if (or *both-regimes?* (not= (law-kind law) :analytical))
    law
    (if (< (rand) 0.06)
      (assoc law :domain (rand-nth [:bound :unbound]))
      law)))

(declare fast-mutate-law fast-crossover-single-laws)

(defn- mutate-analytical-gp [ind]
  (if (and *guess-mutations?* (< (rand) 0.45))
    (guess-mutate-individual ind)
    (random-block-mutate-individual ind)))

(defn gp-mutate-analytical
  "One GP edit on an analytical law (legacy map); no restart or MCTS."
  [ind]
  (validate-and-repair-law
   (if (template-unbound-arms-active?)
     (pair-mutate-analytical-individual ind)
     (mutate-analytical-gp ind))))

(defn- mcts-mutate-roll? []
  (and *mcts-mutate?*
       (seq *mcts-mutate-datasets*)
       (pos? (long *mcts-mutate-simulations*))
       (< (rand) (double *mcts-mutate-rate*))))

(defn- try-mcts-mutate-analytical [ind]
  (try
    (let [search-fn (requiring-resolve 'evophy.mcts/search-repair-individual)
          stats-fn  (requiring-resolve 'evophy.mcts/last-search-result)
          _repaired (search-fn ind (vec *mcts-mutate-datasets*)
                                *mcts-mutate-simulations*
                                :quiet? true)]
      (when (:improved? (stats-fn))
        (:individual (stats-fn))))
    (catch Throwable _ nil)))

(defn- mutate-analytical-individual [ind]
  (let [ind (template-unbound-arms-law ind)]
    (cond
      (mcts-mutate-roll?)
      (or (try-mcts-mutate-analytical ind)
          (if (template-unbound-arms-active?)
            (pair-mutate-analytical-individual ind)
            (mutate-analytical-gp ind)))

      (template-unbound-arms-active?)
      (pair-mutate-analytical-individual ind)

      :else (mutate-analytical-gp ind))))

(defn- mutate-law [law]
  (mutate-analytical-domain
   (law-from-legacy
    (cond
      *fast-breeding?* (fast-mutate-law (law->legacy law))

      (= :analytical (law-kind law))
      (mutate-analytical-individual (law->legacy law))

      (and *guess-mutations?* (< (rand) 0.45))
      (guess-mutate-individual (law->legacy law))

      :else (random-block-mutate-individual (law->legacy law))))))

(defn mutate-individual [ind]
  (let [restart-p (if *stagnation-escape?* 0.4 0.2)]
    (-> (cond
          (< (rand) restart-p) (random-individual)
          (seq (:laws ind))
          (let [i (rand-int (count (:laws ind)))]
            (assoc ind :laws (assoc (vec (:laws ind)) i
                                    (mutate-law (nth (:laws ind) i)))))
          *fast-breeding?* (individual-with-law ind (fast-mutate-law (first-law-legacy ind)))
          (= :analytical (primary-strategy-label ind))
          (individual-with-law ind (mutate-analytical-individual (first-law-legacy ind)))
          (and *guess-mutations?* (< (rand) 0.45)) (guess-mutate-individual ind)
          :else (random-block-mutate-individual ind))
        validate-and-repair)))

;; ── Cross-pollination (GP crossover) ─────────────────────────────────────────

(defn- all-subtrees
  "Flat list of every node in an expression tree (the tree itself plus every subtree)."
  [expr]
  (if (coll? expr)
    (cons expr (mapcat all-subtrees (rest expr)))
    [expr]))

(defn- splice-subtree
  "Walk expr, replacing a random node with donor-sub.
   p = base probability of replacement at the current node; doubles at each level
   so the walk always terminates.  Returns a normalised expression."
  [expr donor-sub p]
  (if (< (rand) (min 1.0 p))
    donor-sub
    (if (coll? expr)
      (let [op   (first expr)
            args (mapv #(splice-subtree % donor-sub (* p 2.0)) (rest expr))]
        (normalize-expr (apply list op args)))
      expr)))

(defn crossover-expr
  "GP subtree crossover: graft a random subtree of donor into a random position in recipient."
  [recipient donor]
  (if (or (expr-too-large? recipient) (expr-too-large? donor))
    (normalize-expr recipient)
    (let [donor-nodes (filterv #(or (coll? %) (symbol? %)) (all-subtrees donor))
          donor-sub   (rand-nth (if (seq donor-nodes) donor-nodes [donor]))]
      (normalize-expr (splice-subtree recipient donor-sub 0.25)))))

(defn- fast-crossover-single-laws
  "Light crossover during escape: per-slot parent pick, no block sync."
  [ind-a ind-b]
  (when (= (law-kind ind-a) (law-kind ind-b))
    (let [strategy (law-kind ind-a)
          ks       (individual-expr-keys ind-a)
          domain   (when (and (= strategy :analytical) (not *both-regimes?*))
                     (let [pick (analytical-law-domain (if (< (rand) 0.5) ind-a ind-b))]
                       (if (= pick :any) (rand-nth [:bound :unbound]) pick)))]
      (cond-> (into {:strategy strategy}
                    (map (fn [k] [k (if (< (rand) 0.5) (get ind-a k) (get ind-b k))]) ks))
        domain (assoc :domain domain)))))

(defn- fast-mutate-law [law]
  (let [leg      (if (:strategy law) law (law->legacy law))
        strategy (law-kind leg)
        [ks vars] (case strategy
                    :analytical   [analytical-expr-keys analytical-vars]
                    :differential [differential-expr-keys differential-vars]
                    :conserved    [[conserved-expr-key] conserved-vars]
                    [nil nil])]
    (if (seq ks)
      (let [k (rand-nth ks)]
        (update leg k #(mutate % vars)))
      leg)))

(defn- crossover-single-laws
  "Cross-pollinate two same-strategy individuals via functional blocks.
   Each expression slot is independently either:
     - kept from parent A           (20%)
     - taken whole from parent B    (20%)
     - block crossover A ← B        (45%)  — swap/append blocks from donor
     - subtree splice (legacy GP)   (15%)
   Returns nil if strategies differ."
  [ind-a ind-b]
  (when (= (:strategy ind-a) (:strategy ind-b))
    (let [strategy (:strategy ind-a)
          ind-a    (sync-block-genome ind-a)
          ind-b    (sync-block-genome ind-b)
          ks       (individual-expr-keys ind-a)
          vars     (case strategy
                     :analytical   analytical-vars
                     :differential differential-vars
                     :conserved    conserved-vars
                     nil)]
      (when (seq ks)
        (letfn [(cross-slot [k]
                  (let [ba (get-in ind-a [:expr-blocks k])
                        bb (get-in ind-b [:expr-blocks k])
                        r  (rand)]
                    (cond
                      (< r 0.20) ba
                      (< r 0.40) bb
                      (< r 0.85) (block-crossover-blocks ba bb)
                      :else (decompose-expr-to-blocks
                              (crossover-expr (get ind-a k) (get ind-b k))))))]
          (let [child-blocks (into {} (map (fn [k] [k (cross-slot k)]) ks))
                domain     (when (and (= strategy :analytical) (not *both-regimes?*))
                             (let [pick (analytical-law-domain (if (< (rand) 0.5) ind-a ind-b))]
                               (if (= pick :any) (rand-nth [:bound :unbound]) pick)))
                child        (cond-> {:strategy strategy :expr-blocks child-blocks}
                               domain (assoc :domain domain))
                child        (-> child (sync-block-genome :from-blocks? true))
                mutated      (reduce (fn [acc k]
                                       (update acc k #(mutate-expr-via-blocks % vars k)))
                                     child ks)]
            (validate-and-repair-law mutated)))))))

(defn crossover-individuals
  "Cross-pollinate two individuals. Composite genomes cross law-by-law when kinds match."
  [ind-a ind-b]
  (if *fast-breeding?*
    (cond
      (and (seq (:laws ind-a)) (seq (:laws ind-b)))
      (when (= (individual-law-kinds ind-a) (individual-law-kinds ind-b))
        (let [laws (filterv some?
                          (mapv (fn [la lb]
                                  (when-let [c (fast-crossover-single-laws (law->legacy la)
                                                                            (law->legacy lb))]
                                    (law-from-legacy c)))
                                (:laws ind-a) (:laws ind-b)))]
          (when (seq laws) (laws-individual laws))))

      :else
      (when-let [child (fast-crossover-single-laws (first-law-legacy ind-a)
                                                   (first-law-legacy ind-b))]
        (individual-with-law ind-a child)))
    (cond
      (and (seq (:laws ind-a)) (seq (:laws ind-b)))
      (when (= (individual-law-kinds ind-a) (individual-law-kinds ind-b))
        (laws-individual
         (mapv (fn [la lb]
                 (law-from-legacy (crossover-single-laws (law->legacy la) (law->legacy lb))))
               (:laws ind-a) (:laws ind-b))))

      :else
      (crossover-single-laws ind-a ind-b))))


(def default-population-file "data/population.edn")
(def checkpoint-version 10)
(def ^:private history-version 1)

(defn history-path-for [population-path]
  (clojure.string/replace population-path #"\.edn$" "-history.edn"))

(defn save-history!
  "Write per-generation fitness stats alongside the population checkpoint."
  [path history]
  (let [f (java.io.File. path)]
    (.mkdirs (.getParentFile f))
    (spit path (pr-str {:version history-version :history (vec history)}))))

(defn load-history
  "Load fitness history from path, or nil if missing / incompatible version."
  [path]
  (let [f (java.io.File. path)]
    (when (.exists f)
      (try
        (let [{:keys [version history]} (edn/read-string (slurp path))]
          (when (= version history-version) (vec history)))
        (catch Exception _ nil)))))

(defn- population-fitness-stats
  "Stats computed from individuals that already carry a :fitness value (the elite tier)."
  [population]
  (let [fits (vec (sort (keep :fitness population)))]
    (if (empty? fits)
      {:best 0.0 :mean 0.0 :median 0.0 :worst 0.0 :n 0}
      (let [n (count fits)
            mid (quot n 2)
            median (if (odd? n)
                     (nth fits mid)
                     (/ (+ (nth fits (dec mid)) (nth fits mid)) 2.0))]
        {:best   (last fits)
         :mean   (/ (apply + fits) (double n))
         :median (double median)
         :worst  (first fits)
         :n      n}))))

(defn- law-for-save [law]
  (let [leg  (sync-block-genome (law->legacy law))
        norm (fn [k] [k (normalize-expr (get leg k))])]
    (into (cond-> {:kind (law-kind leg) :expr-blocks (:expr-blocks leg)}
            (coords/law-chart leg) (assoc :chart (coords/law-chart leg))
            (and (= (law-kind leg) :analytical)
                 (not= (analytical-law-domain leg) :any))
            (assoc :domain (analytical-law-domain leg)))
          (map norm (individual-expr-keys leg)))))

(defn- individual-for-save [ind]
  {:laws (mapv law-for-save (individual-laws ind))})

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
          (when (contains? #{7 8 9 10} version)
            {:population (mapv (fn [ind]
                                 (-> (if (:laws ind)
                                       ind
                                       (migrate-individual-to-laws ind))
                                     sync-individual))
                               population)
             :generations-run (long (or generations-run 0))
             :population-size population-size}))
        (catch Exception e
          (timbre/warn "Could not load population from" path "-" (.getMessage e))
          nil)))))

(defn normalize-population-size
  "Sanitize, pad with random individuals, or trim so size matches the configured target.
   Also runs validate-and-repair on every loaded individual so degenerate genomes
   (e.g. (- py py) = 0 from old checkpoints) are cleaned up before evolution starts."
  [population target-size]
  (let [sanitized (mapv validate-and-repair population)
        n         (count sanitized)]
    (cond
      (zero? target-size) []
      (<= n target-size) (into sanitized (repeatedly (- target-size n) random-individual))
      :else (subvec sanitized 0 target-size))))

(defn- splice-seeds
  "Replace the last (count seeds) slots in population with physics-seed individuals,
   so seeds always enter the first generation regardless of checkpoint state."
  [population seeds]
  (let [pop   (vec population)
        n     (count pop)
        ns    (min (count seeds) n)
        seeds (vec (take ns seeds))]
    (into (subvec pop 0 (- n ns)) seeds)))

(defn resolve-initial-population
  [{:keys [fresh? seed? path population-size strategy de-driven? both-regimes?]}]
  (let [seeds (when seed?
                (cond
                  (and de-driven? (= strategy :analytical) both-regimes?)
                  (both-regimes-analytical-seeds)

                  (and de-driven? (= strategy :analytical))
                  (filterv #(= (:strategy %) :analytical) physics-seeds)

                  de-driven?
                  de-driven-composite-seeds

                  :else
                  (cond->> physics-seeds
                    strategy (filterv #(= (:strategy %) strategy)))))
        base  (if fresh?
                {:population (vec (repeatedly population-size random-individual))
                 :generations-run 0
                 :resumed? false}
                (if-let [{:keys [population generations-run]} (load-population path)]
                  {:population (vec (keep identity
                                          (map #(or (strip-differential-laws %)
                                                    (when de-driven? (random-individual)))
                                               population)))
                   :generations-run generations-run
                   :resumed? true}
                  {:population (vec (repeatedly population-size random-individual))
                   :generations-run 0
                   :resumed? false}))]
    (if (seq seeds)
      (do (println (str "Injecting " (count seeds) " physics seed bundle(s) into initial population."))
          (update base :population splice-seeds seeds))
      base)))

(def default-mcts-simulations 64)
(def default-mcts-inject 5)
(def default-mcts-repair-simulations 48)
(def default-mcts-repair-inject 1)

(defn parse-args
  [args]
  (loop [opts         {:fresh? false
               :seed? false
               :path default-population-file
               :generations 50
               :population-size 50
               :mcts-simulations default-mcts-simulations
               :mcts-until-stop false
               :mcts-inject default-mcts-inject
               :mcts-repair? true
               :mcts-repair-simulations default-mcts-repair-simulations
               :mcts-repair-inject default-mcts-repair-inject
               :mcts-mutate? true
               :mcts-mutate-rate default-mcts-mutate-rate
               :mcts-mutate-simulations default-mcts-mutate-simulations
               :template-unbound-arms? true
               :template-conic-unbound? true
               :analytical-blocks? true
               :mcts? true
               :prompt-each-generation false
               :fitness-mode :random
               :scenario-samples 32
               :fitness-aggregate :min
               :fitness-percentile 10
               :scenario-seed nil
               :guess-mutations? true
               :de-driven? false
               :both-regimes? false
               :domain-filter? false
               :primitive-tier :auto
               :strategy nil}        ; nil = all strategies; or :analytical/:differential/:conserved
         xs args]
    (if (empty? xs)
      (assoc opts :fitness-context
             (make-fitness-context :mode (:fitness-mode opts)
                                   :evaluation (if (:de-driven? opts)
                                                 :de-driven
                                                 :data-driven)
                                   :sample-count (:scenario-samples opts)
                                   :aggregate (:fitness-aggregate opts)
                                   :percentile (:fitness-percentile opts)
                                   :seed (:scenario-seed opts)))
      (let [[a & more] xs]
        (case a
          "--fresh" (recur (assoc opts :fresh? true) more)
          "--seed"  (recur (assoc opts :seed? true) more)
          "--no-mcts" (recur (assoc opts :mcts? false) more)
          "--no-mcts-repair" (recur (assoc opts :mcts-repair? false) more)
          "--mcts-repair" (recur (assoc opts :mcts-repair? true) more)
          "--no-mcts-mutate" (recur (assoc opts :mcts-mutate? false) more)
          "--mcts-mutate" (recur (assoc opts :mcts-mutate? true) more)
          "--no-template-unbound" (recur (assoc opts :template-unbound-arms? false) more)
          "--template-unbound" (recur (assoc opts :template-unbound-arms? true) more)
          "--no-template-conic-unbound" (recur (assoc opts :template-conic-unbound? false) more)
          "--template-conic-unbound" (recur (assoc opts :template-conic-unbound? true) more)
          "--no-analytical-blocks" (recur (assoc opts :analytical-blocks? false) more)
          "--analytical-blocks" (recur (assoc opts :analytical-blocks? true) more)
          "--mcts-until-stop" (recur (assoc opts :mcts-until-stop true) more)
          "--prompt-each-generation" (recur (assoc opts :prompt-each-generation true) more)
          "--fixed-scenarios" (recur (assoc opts :fitness-mode :fixed) more)
          "--random-scenarios" (recur (assoc opts :fitness-mode :random) more)
          "--scenario-samples" (recur (assoc opts :scenario-samples (Long/parseLong (first more))) (rest more))
          "--fitness-aggregate" (recur (assoc opts :fitness-aggregate (keyword (first more))) (rest more))
          "--fitness-percentile" (recur (assoc opts :fitness-percentile (Long/parseLong (first more))) (rest more))
          "--scenario-seed" (recur (assoc opts :scenario-seed (Long/parseLong (first more))) (rest more))
          "--de-driven" (recur (assoc opts :de-driven? true) more)
          "--both-regimes" (recur (assoc opts :both-regimes? true) more)
          "--domain-filter" (recur (assoc opts :domain-filter? true) more)
          "--no-guess" (recur (assoc opts :guess-mutations? false) more)
          "--primitive-tier"
          (recur (assoc opts :primitive-tier
                         (let [v (first more)]
                           (if (= v "auto") :auto (Long/parseLong v))))
                 (rest more))
          "--mcts-simulations" (recur (assoc opts :mcts-simulations (Long/parseLong (first more))) (rest more))
          "--mcts-inject" (recur (assoc opts :mcts-inject (Long/parseLong (first more))) (rest more))
          "--mcts-repair-simulations" (recur (assoc opts :mcts-repair-simulations (Long/parseLong (first more))) (rest more))
          "--mcts-repair-inject" (recur (assoc opts :mcts-repair-inject (Long/parseLong (first more))) (rest more))
          "--mcts-mutate-rate" (recur (assoc opts :mcts-mutate-rate (Double/parseDouble (first more))) (rest more))
          "--mcts-mutate-simulations" (recur (assoc opts :mcts-mutate-simulations (Long/parseLong (first more))) (rest more))
          "--population" (recur (assoc opts :path (first more)) (rest more))
          "--generations" (recur (assoc opts :generations (Long/parseLong (first more))) (rest more))
          "--population-size" (recur (assoc opts :population-size (Long/parseLong (first more))) (rest more))
          "--strategy" (recur (assoc opts :strategy (keyword (first more))) (rest more))
          (throw (ex-info "Unknown argument"
                          {:arg a
                           :hint "--fresh --de-driven --both-regimes --domain-filter --strategy analytical ..."})))))))

(defn- distinct-elites
  "Elite tier for one generation. During escape burst, collapse behaviorally identical clones."
  [scored elite-cap probes behavior-diverse? cache]
  (cond
    *fast-breeding?* (take elite-cap scored)
    behavior-diverse? (take-distinct-by-behavior elite-cap scored probes cache)
    :else (take elite-cap (distinct scored))))

(defn- pick-other-elite [parent elites]
  (let [others (vec (remove #(identical? % parent) elites))]
    (when (seq others) (rand-nth others))))

(defn evolve-generation
  [population fitness-ctx population-size generation-index
   & {:keys [extra-immigrants elite-divisor behavior-probes behavior-diverse-elites?
             behavior-cache score-progress? gen-label repair-immigrants]
      :or {extra-immigrants 0 elite-divisor 5 behavior-diverse-elites? false
           score-progress? false repair-immigrants []}}]
  (let [evaluation  (:evaluation fitness-ctx :data-driven)
        de-driven? (= evaluation :de-driven)
        datasets    (if de-driven?
                      (mapv #(subsample-dataset % (or *fitness-max-states* 48))
                            (reference-datasets))
                      (datasets-for-fitness-context fitness-ctx :generation generation-index))
        phase-states (when de-driven?
                       (phase-states-for-fitness-context fitness-ctx
                                                         :generation generation-index))
        fit-opts    (select-keys fitness-ctx [:aggregate :percentile :evaluation])
        cache       (or behavior-cache (behavior-key-cache))
        n-pop       (count population)
        ;; 10% fresh random immigrants each generation (+ extras during stagnation).
        immigrant-n (+ (max 1 (quot population-size 10)) extra-immigrants)
        elite-cap   (max 1 (quot population-size elite-divisor))
        scored      (vec
                     (binding [*preferred-law-chart* (coords/dominant-chart datasets)
                               *fast-differential-fitness?* true
                               *fast-conserved-fitness?* true
                               *fitness-max-states* 64]
                       (map-indexed
                        (fn [i ind]
                          (when score-progress?
                            (do (println (format "    %s scoring %d/%d (%s)..."
                                                 (or gen-label "gen") (inc i) n-pop
                                                 (name (primary-strategy-label ind))))
                                (flush)))
                          (assoc ind :fitness
                                 (or (call-with-timeout
                                      (or *fitness-timeout-ms* 0)
                                      (name (primary-strategy-label ind))
                                      #(calculate-fitness-scenarios
                                        ind datasets
                                        :aggregate (:aggregate fit-opts)
                                        :percentile (:percentile fit-opts)
                                        :evaluation evaluation
                                        :phase-states phase-states))
                                     0.0)))
                        population)))
        _ (when score-progress?
            (do (println (format "    %s scored; breeding..." (or gen-label "gen")))
                (flush)))
        scored      (sort-by :fitness #(compare %2 %1) scored)
        unique-elites (do
                        (when score-progress?
                          (println (format "    %s elites (%d)..." (or gen-label "gen")
                                           elite-cap))
                          (flush))
                        (distinct-elites scored elite-cap behavior-probes
                                         behavior-diverse-elites? cache))
        n-elites      (count unique-elites)
        breed-slots   (max 0 (- population-size immigrant-n))
        branch        (long (Math/ceil (/ breed-slots (double n-elites))))
        bred          (do
                        (when score-progress?
                          (println (format "    %s mutating %d slots..."
                                           (or gen-label "gen") breed-slots))
                          (flush))
                        (binding [*mcts-mutate-datasets* (when *mcts-mutate?* datasets)]
                          (take breed-slots
                                (mapcat (fn [parent]
                                          (cons parent
                                                (take (dec branch)
                                                      (repeatedly
                                                       #(if-let [other (when (and (> (count unique-elites) 1)
                                                                                   (< (rand) 0.7))
                                                                              (pick-other-elite parent unique-elites))]
                                                          (or (crossover-individuals parent other)
                                                              (mutate-individual parent))
                                                          (mutate-individual parent))))))
                                        unique-elites))))
        repairs     (vec (take immigrant-n repair-immigrants))
        random-n    (max 0 (- immigrant-n (count repairs)))]
    (vec (take population-size
               (concat
                bred
                repairs
                (do
                  (when (or score-progress? (seq repairs))
                    (println (format "    %s immigrants (%d%s)..."
                                     (or gen-label "gen") immigrant-n
                                     (if (seq repairs)
                                       (str ", " (count repairs) " MCTS repair")
                                       "")))
                    (flush))
                  (binding [*fast-immigrants?* (or *fast-immigrants?* (> extra-immigrants 5))
                            *preferred-law-chart* (coords/dominant-chart datasets)]
                    (repeatedly random-n random-individual))))))))

(defn- prompt-continue-evolution? []
  (print "  Enter = next generation, q = stop and save: ")
  (flush)
  (not= "q" (clojure.string/trim (or (read-line) ""))))

(defn -main [& args]
  (timbre/merge-config! {:min-level :warn})
  (let [{:keys [fresh? seed? path generations population-size prompt-each-generation
                fitness-context fitness-mode scenario-samples fitness-aggregate fitness-percentile
                strategy guess-mutations? both-regimes? domain-filter?
                mcts-repair? mcts-repair-simulations mcts-repair-inject
                mcts-mutate? mcts-mutate-rate mcts-mutate-simulations
                template-unbound-arms?
                template-conic-unbound?
                analytical-blocks?
                primitive-tier]}
        (parse-args args)
        de-driven? (= :de-driven (:evaluation fitness-context))
        both-regimes-on? (or both-regimes?
                             (and de-driven? (= strategy :analytical) (not domain-filter?)))
        mcts-repair-on? (and mcts-repair? de-driven? (= strategy :analytical))
        mcts-mutate-on? (and mcts-mutate? de-driven? (= strategy :analytical))
        template-unbound-on? (and template-unbound-arms? both-regimes-on?
                                  de-driven? (= strategy :analytical))
        analytical-blocks-on? (and analytical-blocks? both-regimes-on?
                                   de-driven? (= strategy :analytical))
        template-conic-on? (and template-conic-unbound? template-unbound-on?
                                analytical-blocks-on? (blocks/kepler-conic-valid?))
        report-scenarios default-scenarios
        ref-datasets (scenarios->datasets report-scenarios)
        preferred-chart (coords/dominant-chart ref-datasets)
        primitive-tier-on? (and de-driven? (= strategy :analytical))
        initial-primitive-tier (if primitive-tier-on?
                                 (if (= primitive-tier :auto)
                                   0
                                   (long primitive-tier))
                                 0)
        primitive-tier-atom (atom initial-primitive-tier)]
  (binding [*strategy-filter* strategy
            *guess-mutations?* (if (false? guess-mutations?) false *guess-mutations?*)
            *de-driven-search?* de-driven?
            *both-regimes?* both-regimes-on?
            *template-unbound-arms?* template-unbound-on?
            *template-conic-unbound?* template-conic-on?
            *analytical-blocks?* analytical-blocks-on?
            *primitive-tier* initial-primitive-tier
            *mcts-mutate?* mcts-mutate-on?
            *mcts-mutate-rate* mcts-mutate-rate
            *mcts-mutate-simulations* mcts-mutate-simulations
            *preferred-law-chart* preferred-chart]
  (when (and de-driven? (= strategy :differential))
    (println "warning: --strategy differential with --de-driven scores 0 (DE is already known)"))
  (when de-driven?
    (println
     (if (= strategy :analytical)
       (str "Each individual: analytical laws only (20% + full-orbit fit)"
            (when both-regimes-on?
              (str "; both regimes (per-regime arm fitness, strict e/if on energy)"
                   (when template-unbound-on?
                     (str ", "
                          (if template-conic-on?
                            "hyperbola-conic unbound template"
                            "Taylor unbound template")
                          " + q/p pair mutations"))))
       (if (= strategy :conserved)
         "Each individual: analytical laws + conserved laws (composite; min fitness)"
         "Each individual: analytical laws (motion DE) + conserved laws (invariants along orbits, not DE solutions)"))))
  (when (and mcts-repair-on? de-driven? (= strategy :analytical))
    (println (str "MCTS adversarial repair: up to " mcts-repair-inject
                  " immigrant(s)/gen, " mcts-repair-simulations " sims each")))
  (when (and mcts-mutate-on? de-driven? (= strategy :analytical))
    (println (str "MCTS mutation: " (int (* 100 mcts-mutate-rate))
                  "% of analytical mutates, " mcts-mutate-simulations " sims each")))
  (when template-unbound-on?
    (println (if template-conic-on?
               "Template unbound arms: hyperbola-conic (sinh/cosh F) locked; bound arms mutate in (qx,px)/(qy,py) pairs"
               "Template unbound arms: Taylor locked; bound arms mutate in (qx,px)/(qy,py) pairs")))
  (when analytical-blocks-on?
    (println (str "Analytical blocks (Emmy-validated): "
                  (blocks/catalog-block-count)
                  " catalog laws + kepler-conic inject"
                  " — circle, ellipse, harmonic, conic (sinh/cosh)")))
  (when (and primitive-tier-on? (pos? initial-primitive-tier))
    (println (str "Graded primitives: tier " initial-primitive-tier
                  " — " (clojure.string/join ", " (map name (prims/unlocked-primitive-ops initial-primitive-tier))))))
  (when (and primitive-tier-on? (= primitive-tier :auto))
    (println "Graded primitives: tier auto (unlocks with hall-of-fame fitness)"))
  (let [{:keys [population generations-run resumed?]}
        (resolve-initial-population {:fresh? fresh?
                                     :seed?  seed?
                                     :path path
                                     :population-size population-size
                                     :strategy strategy
                                     :de-driven? de-driven?
                                     :both-regimes? both-regimes-on?})
        initial (normalize-population-size population population-size)
        fit-opts (select-keys fitness-context [:aggregate :percentile :evaluation])
        eval-phase-states (when de-driven?
                            (phase-states-for-fitness-context
                             (assoc fitness-context :seed (or (:seed fitness-context) 42))
                             :generation 0))
        stopped-early? (atom false)
        checkpoint (atom {:pop initial :generations-run generations-run})
        history-path (history-path-for path)
        history (atom (if fresh? [] (or (load-history history-path) [])))
        save-checkpoint!
        (fn [pop gens]
          (reset! checkpoint {:pop pop :generations-run gens})
          (save-population! path pop
                            :generations-run gens
                            :population-size population-size))
        _ (do (save-checkpoint! initial generations-run)
              (println "checkpoint:" path
                       (if resumed? "(resuming)" "(new checkpoint file)")))
        ;; Fixed reference datasets used only for stable per-generation eval — never for selection.
        ;; Hall-of-fame: best individual ever seen by eval fitness. Injected every generation
        ;; so it can never be lost to random-batch variance.
        hall-of-fame  (atom nil)
        ;; Stagnation tracking: count consecutive gens without HoF improvement.
        stagnation    (atom 0)
        escape-stagnation (atom 0) ; flat HoF while in [escape] burst mode
        stagnation-threshold 20   ; gens of flat HoF before escape mode
        stagnation-burst-threshold 50 ; gens before aggressive pool shake-up
        escape-deep-threshold 25  ; [escape] gens without HoF gain → near-full reseed
        hof-target 0.99
        behavior-probes-lite (build-behavior-probes-lite ref-datasets)
        final-pop (reduce
                   (fn [pop gen-idx]
                     (let [prev-hof-fitness (or (:eval-fitness @hall-of-fame) 0.0)
                           stagnating?      (>= @stagnation stagnation-threshold)
                           suboptimal-hof?  (< prev-hof-fitness hof-target)
                           escape?          (and stagnating? suboptimal-hof?)
                           escape-burst?    (and escape? (>= @stagnation stagnation-burst-threshold))
                           escape-deep?     (and escape-burst? (>= @escape-stagnation escape-deep-threshold))
                           _ (when escape?
                               (do (println (format "  gen %d/%d: evolving%s..."
                                                    (inc gen-idx) generations
                                                    (cond escape-deep? " [escape+]"
                                                          escape-burst? " [escape]"
                                                          :else " [stag]")))
                                   (flush)))
                           behavior-cache   (behavior-key-cache)
                           extra-imm        (cond
                                              escape-deep? (quot population-size 4)
                                              escape-burst? (quot population-size 4)
                                              escape? (quot population-size 4)
                                              stagnating? (quot population-size 10)
                                              :else 0)
                           elite-divisor    (cond
                                              escape-deep? 20
                                              escape-burst? 10
                                              :else 5)
                           hof-gk           (when (and escape-burst? @hall-of-fame)
                                              (individual-genome-key (:ind @hall-of-fame)))
                           _ (when escape-deep?
                               (println (if (and (= strategy :analytical) (analytical-blocks-active?))
                                          (format "  [escape+] catalog reseed (%d slots, HoF held back)"
                                                  population-size)
                                          (format "  [escape+] full reseed (%d slots, HoF held back)"
                                                  population-size))))
                           pop-diverse      (binding [*fast-immigrants?* (or escape-burst? escape-deep?)]
                                              (cond
                                                escape-deep?
                                                (escape-deep-reseed-population population-size)

                                                escape-burst?
                                                (let [without-basin (if hof-gk
                                                                      (vec (remove #(= hof-gk
                                                                                       (individual-genome-key %))
                                                                                   pop))
                                                                      pop)
                                                      n-keep (max 0 (quot population-size 5))
                                                      kept   (vec (take n-keep without-basin))
                                                      n-fresh (- population-size (count kept))]
                                                  (into kept (repeatedly n-fresh random-individual)))

                                                :else pop))
                           ;; HoF injected after escape shake-ups — except [escape+] catalog reseed (explore past HoF basin).
                           pop-for-breed    (let [base (let [v (vec pop-diverse)]
                                                         (if (< (count v) population-size)
                                                           (into v (repeatedly (- population-size (count v))
                                                                               random-individual))
                                                           (vec (take population-size v))))]
                                              (if (and @hall-of-fame (not escape-deep?))
                                                (assoc base (dec population-size) (:ind @hall-of-fame))
                                                base))
                           repair-n (when mcts-repair-on?
                                      (cond
                                        escape-deep? 0
                                        escape-burst? (max mcts-repair-inject 3)
                                        escape? (max mcts-repair-inject 2)
                                        stagnating? mcts-repair-inject
                                        :else mcts-repair-inject))
                           repair-source (when (and repair-n @hall-of-fame)
                                           (:ind @hall-of-fame))
                           repair-immigrants (when (and repair-n repair-source)
                                               (try
                                                 ((requiring-resolve 'evophy.mcts/generation-repair-inject)
                                                  repair-source ref-datasets mcts-repair-simulations repair-n)
                                                 (catch Throwable e
                                                   (println "  warning: MCTS repair failed:" (.getMessage e))
                                                   [])))
                           pop' (binding [*stagnation-escape?* (or escape-burst? escape-deep?)
                                          *fast-breeding?* (or escape-burst? escape-deep?)
                                          *fast-immigrants?* (or escape-burst? escape-deep?)
                                          *primitive-tier* @primitive-tier-atom
                                          *fitness-timeout-ms* (when (or escape-burst? escape-deep?)
                                                                 90000)]
                                  (evolve-generation pop-for-breed fitness-context population-size gen-idx
                                                     :extra-immigrants extra-imm
                                                     :elite-divisor elite-divisor
                                                     :behavior-probes behavior-probes-lite
                                                     :behavior-diverse-elites? false
                                                     :behavior-cache behavior-cache
                                                     :score-progress? (or escape-burst? escape-deep?)
                                                     :repair-immigrants (or repair-immigrants [])
                                                     :gen-label (format "gen %d" (inc gen-idx))))
                           gens (+ generations-run (inc gen-idx))
                           {:keys [best mean median]} (population-fitness-stats pop')
                           ;; Re-eval top training elites on fixed reference scenarios (not whole pop).
                           elite-inds (->> pop' (sort-by :fitness #(compare %2 %1)) (take 5))
                           eval-timeout (when (or escape-burst? escape-deep?) 90000)
                           eval-pairs (keep (fn [ind]
                                              (let [s (call-with-timeout
                                                       eval-timeout "eval"
                                                       #(calculate-fitness-scenarios
                                                         ind ref-datasets
                                                         :aggregate (:aggregate fit-opts)
                                                         :percentile (:percentile fit-opts)
                                                         :evaluation (:evaluation fit-opts)
                                                         :phase-states eval-phase-states))]
                                                (when (pos? s) [ind s])))
                                            elite-inds)
                           [best-eval-ind eval-best] (when (seq eval-pairs)
                                                       (apply max-key second eval-pairs))
                           ;; Update hall of fame if this generation produced a new best.
                           _ (when (and eval-best
                                        (or (nil? @hall-of-fame)
                                            (> eval-best (:eval-fitness @hall-of-fame))))
                               (reset! hall-of-fame {:ind best-eval-ind :eval-fitness eval-best}))
                           _ (when (and primitive-tier-on? (= primitive-tier :auto))
                               (reset! primitive-tier-atom
                                       (prims/primitive-tier-for-hof
                                        (or (:eval-fitness @hall-of-fame) eval-best 0.0))))
                           ;; Update stagnation counter.
                           _ (if (> (or (:eval-fitness @hall-of-fame) 0.0) prev-hof-fitness)
                               (do (reset! stagnation 0)
                                   (reset! escape-stagnation 0))
                               (do (swap! stagnation inc)
                                   (when escape-burst?
                                     (swap! escape-stagnation inc))))]
                       (save-checkpoint! pop' gens)
                       (swap! history conj {:gen (inc gen-idx) :total-gen gens
                                            :best best :mean mean :median median
                                            :eval-best eval-best
                                            :hof-best (:eval-fitness @hall-of-fame)})
                       (save-history! history-path @history)
                       (println (format "  gen %d/%d  train-best=%.4f  eval-best=%.4f  hof=%.4f  mean=%.4f%s  saved"
                                        (inc gen-idx) generations best
                                        (or eval-best 0.0)
                                        (or (:eval-fitness @hall-of-fame) 0.0)
                                        mean
                                        (str (when escape-deep? "  [escape+]")
                                             (when (and escape-burst? (not escape-deep?)) "  [escape]")
                                             (when (and escape? (not escape-burst?))
                                               (str "  [stag=" @stagnation "]"))
                                             (when (and stagnating? (not escape?))
                                               (str "  [stag=" @stagnation "]")))))
                       (if (and prompt-each-generation
                                (< gen-idx (dec generations))
                                (not (prompt-continue-evolution?)))
                         (do (reset! stopped-early? true) (reduced pop'))
                         pop')))
                   initial
                   (range generations))
        total-generations (:generations-run @checkpoint)
        final-datasets (if de-driven?
                         ref-datasets
                         (datasets-for-fitness-context fitness-context :generation generations))
        final-probes (build-behavior-probes final-datasets)
        ranked (->> final-pop
                    (map (fn [ind]
                           (assoc ind :fitness (calculate-fitness-scenarios
                                                ind final-datasets
                                                :aggregate (:aggregate fit-opts)
                                                :percentile (:percentile fit-opts)
                                                :evaluation (:evaluation fit-opts)
                                                :phase-states eval-phase-states))))
                    (sort-by :fitness #(compare %2 %1))
                    vec)
        best (first ranked)
        top5 (take-distinct-by-behavior 5 ranked final-probes)]
    (save-checkpoint! final-pop total-generations)
    (println (if resumed? "resumed from" "finished; checkpoint") path)
    (println "total generations (cumulative):" total-generations)
    (println (str "fitness mode: " (name fitness-mode)
                  (when de-driven?
                    " | evaluation: DE-driven (analytical 20%+full-orbit + conserved invariance)")
                  (when-not de-driven?
                    (str " | scenarios/gen: " (if (= fitness-mode :fixed)
                                                (count default-scenarios)
                                                scenario-samples)))
                  " | aggregate: " (name fitness-aggregate)
                  (when (= fitness-aggregate :percentile)
                    (str " | p" fitness-percentile))
                  (when de-driven?
                    (str " | phase-samples/gen: " (:phase-samples fitness-context)))))
    (when (and (not de-driven?) (= fitness-mode :random))
      (println "  (new random scenario batch each generation; use --fixed-scenarios for the 5 named orbits)"))
    (when @stopped-early?
      (println "stopped early (q at prompt) — saved best-so-far population"))
    (println "\nTop 5:")
    (doseq [[i ind] (map-indexed vector top5)]
      (let [fitness (:fitness ind)
            label (primary-strategy-label ind)
            n-laws (count (individual-laws ind))]
        (println (str "  #" (inc i) " " (name label)
                      (when (> n-laws 1) (str " (" n-laws " laws)"))
                      " fitness=" (format "%.4f" fitness)))
        (doseq [[li law] (map-indexed vector (individual-laws ind))]
          (when (> n-laws 1)
            (println (str "    law " (inc li) " (" (name (law-kind law)) "):")))
          (case (law-kind law)
            :analytical
            (do
              (when (analytical-branches-on-energy? law)
                (println (if (analytical-strict-energy-branches? law)
                             "    regime:strict (every slot e/if neg? energy)"
                             "    regime:partial (e/if neg? energy in some slots)")))
              (when (and (not both-regimes-on?)
                         (not= (analytical-law-domain law) :any)
                         (not (analytical-branches-on-energy? law)))
                (println (str "    domain:" (name (analytical-law-domain law)))))
              (let [keys (coords/analytical-expr-keys-for-chart (infer-law-chart law))]
                (doseq [k keys]
                  (println (str "    " (name k) ":" (pr-str (normalize-expr (get law k))))))))
            :differential
            (doseq [k differential-expr-keys]
              (println (str "    " (name k) ":" (pr-str (normalize-expr (get law k))))))
            :conserved
            (println (str "    c-expr:" (pr-str (normalize-expr (get law conserved-expr-key)))))
            (println "    (unknown law kind)")))))
    (println "\nBest on fixed reference scenarios (MSE):")
    (doseq [scenario report-scenarios]
      (let [dataset (scenario-data scenario)]
        (doseq [[li law] (map-indexed vector (individual-laws best))]
          (let [metrics (evaluate-predictions (law->legacy law) dataset)
                regime  (or (:regime dataset) (scenario-regime dataset))]
            (println (str "  " (name (:id scenario))
                          " [" (name regime) "]"
                          (when (> (count (individual-laws best)) 1)
                            (str " law-" (inc li)))
                          " " (name (:strategy metrics))
                          (when (= (:strategy metrics) :analytical)
                            (cond
                              (analytical-strict-energy-branches? (law->legacy law))
                              " regime=per-arm"
                              (analytical-branches-on-energy? (law->legacy law))
                              " regime=branched"
                              both-regimes-on?
                              " regime=both"
                              :else
                              (str " domain=" (name (analytical-law-domain (law->legacy law))))))
                          " mse=" (format-mse (:mse metrics)))))))))))))
