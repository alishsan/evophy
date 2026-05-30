
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

;; ── Exact Keplerian solver ────────────────────────────────────────────────
;; Solves Kepler's equation via Newton's method so generated data has zero
;; orbital drift regardless of how many periods are integrated.

(defn- solve-kepler-elliptic
  "Solve M = u − e·sin(u) for eccentric anomaly u (elliptic orbit, 0 ≤ e < 1)."
  ^double [^double e ^double M]
  (let [twopi (* 2.0 Math/PI)
        M     (- M (* twopi (Math/floor (/ M twopi))))]
    (loop [u (double M) n 0]
      (let [du (/ (- u (* e (Math/sin u)) M)
                  (- 1.0 (* e (Math/cos u))))]
        (if (or (< (Math/abs du) 1e-12) (> n 60))
          u
          (recur (- u du) (inc n)))))))

(defn- solve-kepler-hyperbolic
  "Solve Mh = e·sinh(F) − F for hyperbolic anomaly F (hyperbolic orbit, e > 1)."
  ^double [^double e ^double Mh]
  (loop [F (double Mh) n 0]
    (let [dF (/ (- (* e (Math/sinh F)) F Mh)
                (- (* e (Math/cosh F)) 1.0))]
      (if (or (< (Math/abs dF) 1e-12) (> n 60))
        F
        (recur (- F dF) (inc n))))))

(defn- kepler-state-at-t
  "Exact Keplerian (qx, qy, px, py) at time t via orbital elements.
   Handles elliptic (E < 0) and hyperbolic (E > 0).
   Returns nil for degenerate orbits (radial: |L| < 1e-8; near-parabolic: |e−1| < 1e-5)
   so callers can fall back to numerical integration."
  [m alpha q0x q0y p0x p0y t]
  (let [q0x   (double q0x) q0y (double q0y)
        p0x   (double p0x) p0y (double p0y)
        m     (double m)   alpha (double alpha) t (double t)
        r0    (Math/sqrt (+ (* q0x q0x) (* q0y q0y)))
        L     (- (* q0x p0y) (* q0y p0x))         ; scalar angular momentum
        E     (grav2d-energy m alpha {:qx q0x :qy q0y :px p0x :py p0y})
        ;; Laplace-Runge-Lenz vector: e_vec = (p×L)/(mα) − r̂
        ;; In 2D with L scalar: p×L = (py·L, −px·L)
        ex    (- (/ (* p0y L) (* m alpha)) (/ q0x r0))
        ey    (- (/ (* (- p0x) L) (* m alpha)) (/ q0y r0))
        e     (Math/sqrt (+ (* ex ex) (* ey ey)))]
    (when (and (> (Math/abs L) 1e-8)
               (> (Math/abs (- e 1.0)) 1e-5))
      (let [omega  (Math/atan2 ey ex)          ; angle of periapsis from +x axis
            cos-om (Math/cos omega)
            sin-om (Math/sin omega)
            ;; Rotate initial position into orbital frame (periapsis along +x')
            x0p    (+ (* q0x cos-om) (* q0y sin-om))
            y0p    (+ (* (- q0x) sin-om) (* q0y cos-om))]
        (if (neg? E)
          ;; ── Elliptic orbit ──────────────────────────────────────────────
          ;; Parameterisation: x'=a(cos u − e), y'=b·sin u, r=a(1−e·cos u)
          ;; Time equation: M = u − e·sin u,  M = n(t−t₀), n = √(α/ma³)
          (let [a   (/ alpha (* -2.0 E))
                b   (* a (Math/sqrt (max 0.0 (- 1.0 (* e e)))))
                n   (Math/sqrt (/ alpha (* m a a a)))
                u0  (Math/atan2 (/ y0p b) (+ (/ x0p a) e))
                M0  (- u0 (* e (Math/sin u0)))
                u   (solve-kepler-elliptic e (+ M0 (* n t)))
                cu  (Math/cos u) su (Math/sin u)
                xp  (* a (- cu e))
                yp  (* b su)
                vf  (/ (* a n) (- 1.0 (* e cu)))   ; common velocity factor
                vxp (* (- su) vf)
                vyp (* (/ b a) cu vf)
                qx  (- (* xp cos-om) (* yp sin-om))
                qy  (+ (* xp sin-om) (* yp cos-om))
                vx  (- (* vxp cos-om) (* vyp sin-om))
                vy  (+ (* vxp sin-om) (* vyp cos-om))]
            {:qx qx :qy qy :px (* m vx) :py (* m vy)})
          ;; ── Hyperbolic orbit ─────────────────────────────────────────────
          ;; Parameterisation: x'=a(e−cosh F), y'=bh·sinh F, r=a(e·cosh F−1)
          ;; Time equation: Mh = e·sinh F − F,  Mh = n·t,  n = √(α/ma³)
          (let [a   (/ alpha (* 2.0 E))
                bh  (* a (Math/sqrt (max 0.0 (- (* e e) 1.0))))
                n   (Math/sqrt (/ alpha (* m a a a)))
                sh0 (/ y0p bh)
                F0  (Math/log (+ sh0 (Math/sqrt (+ 1.0 (* sh0 sh0)))))  ; asinh(sh0)
                M0h (- (* e (Math/sinh F0)) F0)
                F   (solve-kepler-hyperbolic e (+ M0h (* n t)))
                chF (Math/cosh F) shF (Math/sinh F)
                xp  (* a (- e chF))
                yp  (* bh shF)
                vf  (/ (* a n) (- (* e chF) 1.0))
                vxp (* (- shF) vf)
                vyp (* (/ bh a) chF vf)
                qx  (- (* xp cos-om) (* yp sin-om))
                qy  (+ (* xp sin-om) (* yp cos-om))
                vx  (- (* vxp cos-om) (* vyp sin-om))
                vy  (+ (* vxp sin-om) (* vyp cos-om))]
            {:qx qx :qy qy :px (* m vx) :py (* m vy)}))))))

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
  "Scenario map plus integrated :data trajectory (includes :m :alpha for fitness)."
  [scenario]
  (let [{:keys [m alpha q0x q0y p0x p0y dt steps]} scenario]
    (assoc scenario :data (vec (generate-data m alpha q0x q0y p0x p0y dt steps)))))

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

  :mode — :random (default) or :fixed (uses [[default-scenarios]]).
  :sample-count — random scenarios per evaluation batch (default 32).
  :aggregate — :min (worst scenario) or :percentile.
  :percentile — fraction in [0,1] or whole percent e.g. 10 for 10th percentile (default 0.1).
  :seed — optional; with :generation offset, stabilizes per-generation batches across a run."
  [& {:keys [mode sample-count aggregate percentile seed]
      :or {mode :random
           sample-count 32
           aggregate :min
           percentile 0.1}}]
  (let [pct (if (<= percentile 1) (double percentile) (/ (double percentile) 100.0))]
    {:mode (keyword mode)
     :sample-count (long sample-count)
     :aggregate (keyword aggregate)
     :percentile pct
     :seed (when seed (long seed))
     :datasets (when (= :fixed (keyword mode))
                 (scenarios->datasets default-scenarios))}))

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

(defn- compile-state-fn [expr]
  (binding [*ns* (the-ns 'evophy.core)]
    (eval `(fn [~'t ~'q0x ~'q0y ~'p0x ~'p0y ~'m ~'alpha]
             ~(-> expr simplify-expr rewrite-div-in-expr)))))

(defn- compile-rate-fn [expr]
  (binding [*ns* (the-ns 'evophy.core)]
    (eval `(fn [~'qx ~'qy ~'px ~'py ~'m ~'alpha]
             (let [~'r2 (+ (* ~'qx ~'qx) (* ~'qy ~'qy))
                   ~'r  (Math/sqrt (max ~'r2 1e-24))
                   ~'r3 (* ~'r ~'r2)]
               ~(-> expr simplify-expr rewrite-div-in-expr))))))

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
     horizon-times)))

(defn- predict-step-rates
  "One Euler step from current (qx,qy,px,py) using evolved rate laws."
  [{:keys [dqx dqy dpx dpy]} {:keys [qx qy px py]} dt m alpha]
  (let [md (double m)
        ad (double alpha)
        dqx (safe-double (dqx qx qy px py md ad))
        dqy (safe-double (dqy qx qy px py md ad))
        dpx (safe-double (dpx qx qy px py md ad))
        dpy (safe-double (dpy qx qy px py md ad))
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
            ;; r, r2, r3 implicitly cover both qx and qy (r = sqrt(qx²+qy²))
            (let [position-aliases (when (#{`qx 'qx `qy 'qy} sym) derived-vars)]
              (some (fn [k]
                      (let [expr (get ind k)]
                        (or (expr-uses-symbol? expr sym)
                            (some #(expr-uses-symbol? expr %) position-aliases))))
                    expr-keys)))
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
                [(quant-double (safe-double (qx-fn t q0x q0y p0x p0y m alpha)))
                 (quant-double (safe-double (qy-fn t q0x q0y p0x p0y m alpha)))
                 (quant-double (safe-double (px-fn t q0x q0y p0x p0y m alpha)))
                 (quant-double (safe-double (py-fn t q0x q0y p0x p0y m alpha)))]))]))
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

(defn normalize-expr
  "GP trees must be plain lists for eval/print; mutate used to leave LazySeq subtrees."
  [expr]
  (cond
    (instance? clojure.lang.LazySeq expr) (normalize-expr (doall expr))
    (and (sequential? expr) (not (map? expr)))
    (apply list (map normalize-expr expr))
    :else expr))

(def ^:private latex-symbol-names
  '{t "t" m "m" alpha "\\alpha"
    q0x "q_{0x}" q0y "q_{0y}" p0x "p_{0x}" p0y "p_{0y}"
    qx "q_x" qy "q_y" px "p_x" py "p_y"})

(def ^:private analytical-equation-specs
  [[:qx-expr :qx "q_x(t)" "q_x(t)"]
   [:qy-expr :qy "q_y(t)" "q_y(t)"]
   [:px-expr :px "p_x(t)" "p_x(t)"]
   [:py-expr :py "p_y(t)" "p_y(t)"]])

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
          e/div (let [a (child (first args)) b (child (second args))]
                  (if latex? (str "\\frac{" a "}{" b "}") (str "(" a " / " b ")")))
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
        predicted-next (when dt (predict-step-rates fns s0 dt m alpha))]
    {:index i
     :t (:t s0)
     :state (select-keys s0 [:qx :qy :px :py])
     :rates rates
     :dt dt
     :predicted-next (select-keys predicted-next [:qx :qy :px :py])
     :actual-next (when s1 (select-keys s1 [:qx :qy :px :py]))}))

(defn individual->equations
  "Readable equations for a genome plus optional numeric spot-check.

  individual — map with :strategy and expression keys.
  scenario — nil, scenario :id keyword, or full scenario map (see resolve-scenario).
  opts:
    :format — :plain (default) or :latex
    :sample-t — analytical only: evaluate all four laws at this t
    :sample-index — differential only: trajectory index for state/rates (default: midpoint)

  Returns {:strategy :format :scenario-id :scenario-params :equations :metrics :sample}."
  [individual & {:keys [scenario format sample-t sample-index]
                 :or {format :plain}}]
  (let [sc (resolve-scenario scenario)
        dataset (scenario-data sc)
        specs (case (:strategy individual)
                :analytical analytical-equation-specs
                :differential differential-equation-specs
                nil)
        _ (when (nil? specs)
            (throw (ex-info "Unknown individual strategy" {:strategy (:strategy individual)})))
        equations (build-equation-map individual specs format)
        metrics (evaluate-predictions individual dataset)
        sample (case (:strategy individual)
                 :analytical
                 (when sample-t
                   (sample-analytical-at-t individual dataset sample-t))
                 :differential
                 (sample-differential-at-index individual dataset sample-index))]
    {:strategy (:strategy individual)
     :format format
     :scenario-id (:id sc)
     :scenario-params (select-keys sc [:m :alpha :q0x :q0y :p0x :p0y :dt :steps])
     :equations equations
     :metrics metrics
     :sample sample}))

(def ops '[+ - * e/square e/sin e/cos e/div e/sqrt])
(def analytical-vars (vec (concat '[t] ic-vars param-vars)))
(def differential-vars (vec (concat state-vars param-vars derived-vars)))
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

(defn calculate-analytical-fitness [ind dataset]
  (try
    (when (analytical-genome-valid? ind)
      (let [{:keys [data]} dataset
            D   (char-scale data)
            fns {:qx (compile-state-fn (:qx-expr ind))
                 :qy (compile-state-fn (:qy-expr ind))
                 :px (compile-state-fn (:px-expr ind))
                 :py (compile-state-fn (:py-expr ind))}
            errors (horizon-errors-analytical fns dataset (all-horizon-times data))]
        (fitness-from-errors errors D)))
    (catch Exception _ 0)))

(defn calculate-differential-fitness [ind dataset]
  (try
    (when (genome-valid? ind)
      (let [{:keys [data]} dataset
            D   (char-scale data)
            fns {:dqx (compile-rate-fn (:dqx-expr ind))
                 :dqy (compile-rate-fn (:dqy-expr ind))
                 :dpx (compile-rate-fn (:dpx-expr ind))
                 :dpy (compile-rate-fn (:dpy-expr ind))}
            one-step (one-step-errors-differential fns dataset)
            rollout  (rollout-errors-differential fns dataset)
            errors   {:sq-q (+ (:sq-q one-step) (:sq-q rollout))
                      :sq-p (+ (:sq-p one-step) (:sq-p rollout))
                      :n    (+ (:n one-step)    (:n rollout))}]
        (fitness-from-errors errors D)))
    (catch Exception _ 0)))

(defn calculate-fitness
  "dataset is a scenario map with :data, :m, :alpha."
  [individual dataset]
  (or (case (:strategy individual)
        :analytical (calculate-analytical-fitness individual dataset)
        :differential (calculate-differential-fitness individual dataset)
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

(defn calculate-fitness-scenarios
  "Robust fitness over scenario datasets: :min (worst case) or :percentile (e.g. p10).

  Optional :aggregate and :percentile (see [[aggregate-scenario-fitness]])."
  [individual datasets & {:keys [aggregate percentile] :or {aggregate :min percentile 0.1}}]
  (aggregate-scenario-fitness
   (mapv #(calculate-fitness individual %) datasets)
   :aggregate aggregate
   :percentile percentile))

(defn- expr-subtrees [expr]
  (let [expr (normalize-expr expr)]
    (if (sequential? expr)
      (cons expr (mapcat expr-subtrees (rest expr)))
      [expr])))

(defn- slotted-subtrees [ind]
  (let [keys (case (:strategy ind)
               :analytical analytical-expr-keys
               :differential differential-expr-keys)]
    (for [slot keys, st (expr-subtrees (get ind slot))]
      {:slot slot :motif st})))

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
  (let [donor-nodes (filterv #(or (coll? %) (symbol? %)) (all-subtrees donor))
        donor-sub   (rand-nth (if (seq donor-nodes) donor-nodes [donor]))]
    (normalize-expr (splice-subtree recipient donor-sub 0.25))))

(defn crossover-individuals
  "Cross-pollinate two same-strategy individuals.
   Each expression key is independently either:
     - kept from parent A  (25%)
     - taken whole from parent B  (25%)
     - subtree-spliced A ← random subtree of B  (50%)
   A light mutation pass is applied afterwards.
   Returns nil if strategies differ."
  [ind-a ind-b]
  (when (= (:strategy ind-a) (:strategy ind-b))
    (let [strategy (:strategy ind-a)
          ks       (case strategy
                     :analytical   analytical-expr-keys
                     :differential differential-expr-keys
                     nil)
          vars     (case strategy
                     :analytical   analytical-vars
                     :differential differential-vars
                     nil)]
      (when ks
        (into {:strategy strategy}
              (map (fn [k]
                     (let [ea (get ind-a k) eb (get ind-b k)
                           r  (rand)
                           child (cond
                                   (< r 0.25) ea
                                   (< r 0.50) eb
                                   :else      (crossover-expr ea eb))]
                       [k (mutate child vars)]))
                   ks))))))


(def default-population-file "data/population.edn")
(def checkpoint-version 7)
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
               :mcts-until-stop false
               :mcts-inject default-mcts-inject
               :mcts? true
               :prompt-each-generation false
               :fitness-mode :random
               :scenario-samples 32
               :fitness-aggregate :min
               :fitness-percentile 10
               :scenario-seed nil}
         xs args]
    (if (empty? xs)
      (assoc opts :fitness-context
             (make-fitness-context :mode (:fitness-mode opts)
                                   :sample-count (:scenario-samples opts)
                                   :aggregate (:fitness-aggregate opts)
                                   :percentile (:fitness-percentile opts)
                                   :seed (:scenario-seed opts)))
      (let [[a & more] xs]
        (case a
          "--fresh" (recur (assoc opts :fresh? true) more)
          "--no-mcts" (recur (assoc opts :mcts? false) more)
          "--mcts-until-stop" (recur (assoc opts :mcts-until-stop true) more)
          "--prompt-each-generation" (recur (assoc opts :prompt-each-generation true) more)
          "--fixed-scenarios" (recur (assoc opts :fitness-mode :fixed) more)
          "--random-scenarios" (recur (assoc opts :fitness-mode :random) more)
          "--scenario-samples" (recur (assoc opts :scenario-samples (Long/parseLong (first more))) (rest more))
          "--fitness-aggregate" (recur (assoc opts :fitness-aggregate (keyword (first more))) (rest more))
          "--fitness-percentile" (recur (assoc opts :fitness-percentile (Long/parseLong (first more))) (rest more))
          "--scenario-seed" (recur (assoc opts :scenario-seed (Long/parseLong (first more))) (rest more))
          "--mcts-simulations" (recur (assoc opts :mcts-simulations (Long/parseLong (first more))) (rest more))
          "--mcts-inject" (recur (assoc opts :mcts-inject (Long/parseLong (first more))) (rest more))
          "--population" (recur (assoc opts :path (first more)) (rest more))
          "--generations" (recur (assoc opts :generations (Long/parseLong (first more))) (rest more))
          "--population-size" (recur (assoc opts :population-size (Long/parseLong (first more))) (rest more))
          (throw (ex-info "Unknown argument"
                          {:arg a
                           :hint "--fresh --fixed-scenarios --random-scenarios --scenario-samples N --fitness-aggregate min|percentile --fitness-percentile P --scenario-seed N --no-mcts ..."})))))))

(defn evolve-generation
  [population fitness-ctx population-size generation-index]
  (let [datasets   (datasets-for-fitness-context fitness-ctx :generation generation-index)
        fit-opts   (select-keys fitness-ctx [:aggregate :percentile])
        ;; 10% fresh random immigrants each generation — floor of diversity that
        ;; prevents complete clonal collapse regardless of selection pressure.
        immigrant-n (max 1 (quot population-size 10))
        elite-cap   (max 1 (quot population-size 5))
        scored      (->> population
                         (map (fn [ind]
                                (assoc ind :fitness (calculate-fitness-scenarios
                                                     ind datasets
                                                     :aggregate (:aggregate fit-opts)
                                                     :percentile (:percentile fit-opts)))))
                         (sort-by :fitness #(compare %2 %1)))
        ;; Deduplicate by genome before selecting elites so a clonal population
        ;; can only contribute one elite slot, not all of them.
        unique-elites (take elite-cap (distinct scored))
        n-elites      (count unique-elites)
        breed-slots   (- population-size immigrant-n)
        branch        (long (Math/ceil (/ breed-slots (double n-elites))))]
    (vec (take population-size
               (concat
                ;; Elites + their offspring (crossover or mutation)
                (mapcat (fn [parent]
                          (cons parent
                                (take (dec branch)
                                      (repeatedly
                                       #(if (and (> (count unique-elites) 1) (< (rand) 0.7))
                                          ;; Cross-pollination: breed with a randomly chosen
                                          ;; different elite; fall back to mutation if strategies differ.
                                          (let [other (rand-nth (remove #{parent} unique-elites))]
                                            (or (crossover-individuals parent other)
                                                (mutate-individual parent)))
                                          ;; Pure mutation
                                          (mutate-individual parent))))))
                        unique-elites)
                ;; Fresh random immigrants — keeps the gene pool open
                (repeatedly immigrant-n random-individual))))))

(defn- prompt-continue-evolution? []
  (print "  Enter = next generation, q = stop and save: ")
  (flush)
  (not= "q" (clojure.string/trim (or (read-line) ""))))

(defn -main [& args]
  (timbre/merge-config! {:min-level :warn})
  (let [{:keys [fresh? path generations population-size prompt-each-generation
                fitness-context fitness-mode scenario-samples fitness-aggregate fitness-percentile]}
        (parse-args args)
        report-scenarios default-scenarios
        {:keys [population generations-run resumed?]}
        (resolve-initial-population {:fresh? fresh?
                                     :path path
                                     :population-size population-size})
        initial (normalize-population-size population population-size)
        fit-opts (select-keys fitness-context [:aggregate :percentile])
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
                       (if resumed? "(resuming)" "(new run — use --fresh to ignore existing file)")))
        ;; Fixed reference datasets used only for stable per-generation eval — never for selection.
        ref-datasets (scenarios->datasets default-scenarios)
        ;; Hall-of-fame: best individual ever seen by eval fitness. Injected every generation
        ;; so it can never be lost to random-batch variance.
        hall-of-fame (atom nil)
        final-pop (reduce
                   (fn [pop gen-idx]
                     ;; Inject the hall-of-fame individual into the pool before evolving so it
                     ;; competes for elite slots on the new batch (replacing the last slot).
                     (let [pop-with-hof (if-let [hof @hall-of-fame]
                                          (assoc (vec pop) (dec (count pop)) (:ind hof))
                                          pop)
                           pop' (evolve-generation pop-with-hof fitness-context population-size gen-idx)
                           gens (+ generations-run (inc gen-idx))
                           {:keys [best mean median]} (population-fitness-stats pop')
                           ;; Eval every elite individual on the fixed reference scenarios and
                           ;; pick the best — more robust than trusting training rank alone.
                           elite-inds (filter :fitness pop')
                           eval-pairs (keep (fn [ind]
                                              (let [s (calculate-fitness-scenarios
                                                       ind ref-datasets
                                                       :aggregate (:aggregate fit-opts)
                                                       :percentile (:percentile fit-opts))]
                                                (when (pos? s) [ind s])))
                                            elite-inds)
                           [best-eval-ind eval-best] (when (seq eval-pairs)
                                                       (apply max-key second eval-pairs))
                           ;; Update hall of fame if this generation produced a new best.
                           _ (when (and eval-best
                                        (or (nil? @hall-of-fame)
                                            (> eval-best (:eval-fitness @hall-of-fame))))
                               (reset! hall-of-fame {:ind best-eval-ind :eval-fitness eval-best}))]
                       (save-checkpoint! pop' gens)
                       (swap! history conj {:gen (inc gen-idx) :total-gen gens
                                            :best best :mean mean :median median
                                            :eval-best eval-best
                                            :hof-best (:eval-fitness @hall-of-fame)})
                       (save-history! history-path @history)
                       (println (format "  gen %d/%d  train-best=%.4f  eval-best=%.4f  hof=%.4f  mean=%.4f  saved"
                                        (inc gen-idx) generations best
                                        (or eval-best 0.0)
                                        (or (:eval-fitness @hall-of-fame) 0.0)
                                        mean))
                       (if (and prompt-each-generation
                                (< gen-idx (dec generations))
                                (not (prompt-continue-evolution?)))
                         (do (reset! stopped-early? true) (reduced pop'))
                         pop')))
                   initial
                   (range generations))
        total-generations (:generations-run @checkpoint)
        final-datasets (datasets-for-fitness-context fitness-context :generation generations)
        final-probes (build-behavior-probes final-datasets)
        ranked (->> final-pop
                    (map (fn [ind]
                           (assoc ind :fitness (calculate-fitness-scenarios ind final-datasets
                                                                          :aggregate (:aggregate fit-opts)
                                                                          :percentile (:percentile fit-opts)))))
                    (sort-by :fitness #(compare %2 %1))
                    vec)
        best (first ranked)
        top5 (take-distinct-by-behavior 5 ranked final-probes)]
    (save-checkpoint! final-pop total-generations)
    (println (if resumed? "resumed from" "finished; checkpoint") path)
    (println "total generations (cumulative):" total-generations)
    (println (str "fitness mode: " (name fitness-mode)
                  " | samples/gen: " scenario-samples
                  " | aggregate: " (name fitness-aggregate)
                  (when (= fitness-aggregate :percentile)
                    (str " | p" fitness-percentile))))
    (when (= fitness-mode :random)
      (println "  (new random scenario batch each generation; use --fixed-scenarios for the 5 named orbits)"))
    (when @stopped-early?
      (println "stopped early (q at prompt) — saved best-so-far population"))
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
    (println "\nBest on fixed reference scenarios (MSE):")
    (doseq [scenario report-scenarios]
      (let [dataset (scenario-data scenario)
            metrics (evaluate-predictions best dataset)]
        (println (str "  " (name (:id scenario)) " " (name (:strategy metrics))
                      " mse=" (format "%.6f" (:mse metrics))))))))
