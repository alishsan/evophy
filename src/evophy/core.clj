
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
  :evaluation — :data-driven (trajectory fit) or :de-driven (DE / invariant
  residual in phase space; analytical + conserved only — differential is redundant).
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

(defn- compile-state-fn [expr]
  ;; Derived variables from initial conditions:
  ;;   r0/r02/r03  — initial orbital radius and its powers
  ;;   omega       — circular angular velocity sqrt(α/(m·r₀³)); exact for circular orbit
  ;;   omega-L     — true initial angular velocity L/(m·r₀²) = (q×p)/(m·r₀²); works for all orbits
  ;; Together these let the GP express Taylor corrections, circular rotation, and the
  ;; angular-momentum rotation approximation q(t) ≈ R(Ω_L·t)·q₀.
  (binding [*ns* (the-ns 'evophy.core)]
    (eval `(fn [~'t ~'q0x ~'q0y ~'p0x ~'p0y ~'m ~'alpha]
             (let [~'r02    (+ (* ~'q0x ~'q0x) (* ~'q0y ~'q0y))
                   ~'r0     (Math/sqrt (max ~'r02 1e-24))
                   ~'r03    (* ~'r0 ~'r02)
                   ~'r05    (* ~'r03 ~'r02)
                   ~'r06    (* ~'r03 ~'r03)
                   ~'omega  (Math/sqrt (max 0.0 (/ ~'alpha (* ~'m ~'r03))))
                   ;; Angular velocity from the actual angular momentum L = q×p
                   ~'omega-L (/ (- (* ~'q0x ~'p0y) (* ~'q0y ~'p0x))
                                (* ~'m ~'r02))]
               ~(-> expr simplify-expr rewrite-div-in-expr))))))

(defn- compile-rate-fn [expr]
  (binding [*ns* (the-ns 'evophy.core)]
    (eval `(fn [~'qx ~'qy ~'px ~'py ~'m ~'alpha]
             (let [~'r2 (+ (* ~'qx ~'qx) (* ~'qy ~'qy))
                   ~'r  (Math/sqrt (max ~'r2 1e-24))
                   ~'r3 (* ~'r ~'r2)]
               ~(-> expr simplify-expr rewrite-div-in-expr))))))

(defn- compile-conserved-fn
  "Compile conserved-quantity expression to a fn [qx qy px py m alpha] -> double.
   Derives r, r2, r3 from position — identical binding to compile-rate-fn."
  [expr]
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

(defn- short-horizon-times
  "Times covering only the first `frac` of the trajectory (default 8%).
   Analytical expressions can't predict full Keplerian orbits (transcendental),
   but Taylor expansion approximations work well for short horizons.
   At 8%, the circle arc is ~27° where even 2nd-order Taylor gives <2% error."
  ([data] (short-horizon-times data 0.08))
  ([data frac]
   (let [ts (all-horizon-times data)
         n  (max 1 (int (Math/ceil (* frac (count ts)))))]
     (vec (take n ts)))))

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

(defn expr-uses-t? [expr]
  (expr-uses-symbol? expr 't))

(def analytical-expr-keys    [:qx-expr :qy-expr :px-expr :py-expr])
(def differential-expr-keys  [:dqx-expr :dqy-expr :dpx-expr :dpy-expr])
(def conserved-expr-key      :c-expr)

;; Analytical expressions must reference the initial position (q0x, q0y) so they aren't
;; trivially constant.  We don't mandate p0x/p0y or alpha explicitly: the circular solution
;; is fully determined by position and omega = sqrt(alpha/(m*r03)), so alpha is implicit.
;; Requiring m ensures the GP uses mass (affects the circular frequency).
(def required-analytical-symbols   '[q0x q0y m])
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

(defn- ic-group-variation-ok?
  "True when at least one (m, α) group has ≥2 scenarios whose per-trajectory
   means differ by >5% CoV — the expression must depend on initial conditions,
   not just parameters.  Returns false when no qualifying groups exist."
  [c-fn datasets]
  (let [groups (->> datasets
                    (map (fn [ds] {:m (:m ds) :alpha (:alpha ds) :ds ds}))
                    (group-by (fn [{:keys [m alpha]}] [(double m) (double alpha)]))
                    vals
                    (filter #(>= (count %) 2)))]
    (and (seq groups)
         (some (fn [grp]
                 (let [ms (keep identity (conserved-per-traj-means c-fn (map :ds grp)))
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
         (some #(expr-uses-symbol? sexpr %) state-vars)
         (expr-uses-position? sexpr)
         (some #(expr-uses-symbol? sexpr %) '[px py])
         (try
           (ic-group-variation-ok?
            (compile-conserved-fn expr)
            (reference-datasets))
           (catch Exception _ false)))))

(defn genome-valid?
  "Analytical: each expr uses t; ICs and (m, α) appear across trajectory laws.
   Differential: state coords and (m, α) appear across rate laws.
   Conserved: expression uses at least one dynamic state variable."
  [{:keys [strategy] :as ind}]
  (case strategy
    :analytical   (analytical-genome-valid? ind)
    :differential (symbols-covered-across-exprs? differential-expr-keys ind
                                                 required-differential-symbols)
    :conserved    (conserved-genome-valid? ind)
    false))

(defn individual-genome-key [ind]
  (case (:strategy ind)
    :analytical   (into [:analytical]   (mapv ind analytical-expr-keys))
    :differential (into [:differential] (mapv ind differential-expr-keys))
    :conserved    [:conserved (get ind conserved-expr-key)]
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
    (case (:strategy ind)
      :conserved
      (when (seq conserved)
        (let [c-fn (compile-conserved-fn (:c-expr ind))]
          [:conserved
           (vec (for [{:keys [qx qy px py m alpha]} conserved]
                  (quant-double (safe-double (c-fn qx qy px py m alpha)))))]))
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
      nil)
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

(defn- evaluate-conserved-predictions
  "Metrics for a conserved-quantity individual.
   :mse = max(squared CoV, squared relative RMS of |∇C·f|/|C|) — fails if either
   constancy or the conservation law is violated."
  [ind dataset]
  (try
    (let [c-fn  (compile-conserved-fn (:c-expr ind))
          c-vals (conserved-vals-for-dataset c-fn dataset)
          d-vals (conserved-dot-vals-for-dataset c-fn dataset)
          fvals  (filterv #(Double/isFinite %) c-vals)
          n      (count fvals)]
      (if (< n 3)
        {:strategy :conserved :mse Double/POSITIVE_INFINITY :n-horizons 0}
        (let [mean     (/ (reduce + fvals) n)
              variance (/ (reduce (fn [s v] (let [d (- v mean)] (+ s (* d d)))) 0.0 fvals) n)
              cov-sq   (if (< (Math/abs mean) 1e-8) Double/POSITIVE_INFINITY
                         (/ variance (* mean mean)))
              pairs    (filterv (fn [[c d]] (and (Double/isFinite c) (Double/isFinite d)))
                              (map vector c-vals d-vals))
              dot-sq   (if (empty? pairs) Double/POSITIVE_INFINITY
                         (/ (reduce (fn [s [c d]]
                                      (let [scale (max (Math/abs c) 1e-8)]
                                        (+ s (* (/ d scale) (/ d scale)))))
                                    0.0 pairs)
                            (count pairs)))]
          {:strategy :conserved
           :n-horizons n
           :mse (max cov-sq dot-sq)})))
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
  '{t "t" m "m" alpha "\\alpha"
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
                :analytical   analytical-equation-specs
                :differential differential-equation-specs
                :conserved    conserved-equation-specs
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
                 (sample-differential-at-index individual dataset sample-index)
                 :conserved nil)]
    {:strategy (:strategy individual)
     :format format
     :scenario-id (:id sc)
     :scenario-params (select-keys sc [:m :alpha :q0x :q0y :p0x :p0y :dt :steps])
     :equations equations
     :metrics metrics
     :sample sample}))

(def ops '[+ - * e/square e/sin e/cos e/div e/sqrt])
;; Derived variables pre-computed by compile-state-fn from initial conditions.
;; r0/r02/r03/r05/r06 — initial radius powers (r₀¹,²,³,⁵,⁶)
;; omega             — circular angular velocity sqrt(α/m/r₀³); exact for circular orbit
;; omega-L           — actual initial angular velocity (q×p)/(m·r₀²); valid for any orbit
(def analytical-derived-vars '[r0 r02 r03 r05 r06 omega omega-L])
(def analytical-vars (vec (concat '[t] ic-vars param-vars analytical-derived-vars)))
(def differential-vars (vec (concat state-vars param-vars derived-vars)))

;; Variables for conserved-quantity expressions: current state, params, and
;; derived geometric variables (r, r2, r3 are pre-computed in compile-conserved-fn).
(def conserved-vars (vec (concat state-vars param-vars derived-vars)))

(def constants '[-1.0 -0.5 0.5 1.0 2.0 3.0])


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
         (some #(expr-uses-symbol? sexpr %) state-vars)
         (expr-uses-position? sexpr)
         (some #(expr-uses-symbol? sexpr %) '[px py]))))

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
  [c-fn {:keys [data m alpha]}]
  (let [md (double m) ad (double alpha)]
    (mapv (fn [{:keys [qx qy px py]}]
            (safe-double (c-fn (double qx) (double qy) (double px) (double py) md ad)))
          data)))

(def ^:private conserved-grad-eps 1e-6)

(defn- conserved-dot-at-state
  "Time derivative dC/dt = ∇C·f at one phase-space point, using central differences
   for ∇C and the exact Hamiltonian vector field f = (q̇, ṗ)."
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

(defn- conserved-dot-vals-for-dataset
  "Per-state |dC/dt| = |∇C·f| along a trajectory."
  [c-fn {:keys [data m alpha]}]
  (let [md (double m) ad (double alpha)]
    (mapv (fn [{:keys [qx qy px py]}]
            (let [dot (conserved-dot-at-state c-fn qx qy px py md ad)]
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

(defn- calculate-de-driven-fitness
  "Fitness from the environment DE — not trajectory matching.
   Differential rate laws score 0 (redundant with the known equations of motion)."
  [ind phase-states]
  (case (:strategy ind)
    :analytical   (or (calculate-analytical-equation-fitness ind phase-states) 0.0)
    :conserved    (or (calculate-conserved-equation-fitness ind phase-states) 0.0)
    :differential 0.0
    0.0))

(defn calculate-conserved-fitness
  "Score a conserved-quantity individual on one trajectory dataset.
   Requires low CoV along the orbit AND ∇C·f ≈ 0 (true invariant)."
  [ind dataset]
  (try
    (when (conserved-genome-valid? ind)
      (conserved-trajectory-fitness (compile-conserved-fn (:c-expr ind)) dataset))
    (catch Exception _ 0)))

(def ^:dynamic *fast-conserved-fitness?* false)
(def ^:dynamic *fast-differential-fitness?* false)
(def ^:dynamic *fitness-max-states* nil)
(def ^:dynamic *fitness-timeout-ms* nil)

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
      (let [c-fn  (compile-conserved-fn (:c-expr ind))
            datasets (if-let [ms *fitness-max-states*]
                       (mapv #(subsample-dataset % ms) datasets)
                       datasets)
            per-traj (mapv (fn [ds]
                             (let [vals   (conserved-vals-for-dataset c-fn ds)
                                   finite (filterv #(Double/isFinite %) vals)
                                   n      (count finite)
                                   mean   (when (>= n 1) (/ (reduce + finite) n))]
                               {:vals finite :n n :mean mean :m (:m ds) :alpha (:alpha ds) :ds ds}))
                           datasets)
            traj-means (keep :mean per-traj)
            ic-ok?   (if *fast-conserved-fitness?*
                       true
                       (ic-group-variation-ok? c-fn (reference-datasets)))
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
             (apply min (map #(conserved-trajectory-fitness c-fn (:ds %)) per-traj))))))
    (catch Exception _ 0)))

(def ^:dynamic *strategy-filter*
  "When set to a keyword (:analytical, :differential, :conserved), random-individual
   generates only that strategy.  nil (default) = uniform mix of all three."
  nil)

(def ^:dynamic *fast-immigrants?*
  "When true, conserved immigrants skip expensive IC integration (fitness filters invalid)."
  false)

(def ^:dynamic *de-driven-search?*
  "When true, immigrants and default strategy mix exclude :differential (redundant with known ODE)."
  false)

(defn random-individual []
  (cond
    (and *de-driven-search?* (nil? *strategy-filter*))
    (if (< (rand) 0.5) (random-analytical-individual) (random-conserved-individual))

    *strategy-filter*
    (case *strategy-filter*
      :analytical   (random-analytical-individual)
      :differential (random-differential-individual)
      :conserved    (if *fast-immigrants?*
                      (random-conserved-individual {:strict? false :max-tries 25})
                      (random-conserved-individual)))

    :else
    (let [r (rand)]
      (cond
        (< r 0.34) (random-analytical-individual)
        (< r 0.67) (random-differential-individual)
        :else (if *fast-immigrants?*
                (random-conserved-individual {:strict? false :max-tries 25})
                (random-conserved-individual))))))

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
    :qx-expr '(+ q0x (* (e/div p0x m) t))
    :qy-expr '(+ q0y (* (e/div p0y m) t))
    :px-expr '(+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
    :py-expr '(+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))}

   ;; Second-order: adds ½·F·t²/m correction to positions
   ;; q(t) ≈ q₀ + (p₀/m)t − (α·q₀)/(2m·r₀³)·t²
   {:strategy :analytical
    :qx-expr '(+ q0x (+ (* (e/div p0x m) t) (* (e/div (* (* -0.5 alpha) q0x) (* m r03)) (* t t))))
    :qy-expr '(+ q0y (+ (* (e/div p0y m) t) (* (e/div (* (* -0.5 alpha) q0y) (* m r03)) (* t t))))
    :px-expr '(+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
    :py-expr '(+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))}

   ;; Circular rotation seed (exact for circular orbit, approximation for ellipse/hyperbola):
   ;; q(t) = R(ωt)·q₀  where ω = sqrt(α/(m·r₀³)), pre-computed as `omega`.
   ;; p(t) = m·ω·R(ωt+π/2)·q₀ = m·ω·(−sin,cos)·q₀
   {:strategy :analytical
    :qx-expr '(- (* q0x (e/cos (* omega t))) (* q0y (e/sin (* omega t))))
    :qy-expr '(+ (* q0x (e/sin (* omega t))) (* q0y (e/cos (* omega t))))
    :px-expr '(* (* -1.0 (* m omega)) (+ (* q0x (e/sin (* omega t))) (* q0y (e/cos (* omega t)))))
    :py-expr '(* (* m omega) (- (* q0x (e/cos (* omega t))) (* q0y (e/sin (* omega t)))))}

   ;; Angular-momentum rotation seed: uses the ACTUAL initial angular velocity Ω_L = L/(m·r₀²)
   ;; where L = q₀×p₀ = q0x·p0y − q0y·p0x.  Better than circular-omega for eccentric orbits.
   ;; Momentum follows from d/dt of the rotated position.
   {:strategy :analytical
    :qx-expr '(- (* q0x (e/cos (* omega-L t))) (* q0y (e/sin (* omega-L t))))
    :qy-expr '(+ (* q0x (e/sin (* omega-L t))) (* q0y (e/cos (* omega-L t))))
    :px-expr '(* (* -1.0 (* m omega-L)) (+ (* q0x (e/sin (* omega-L t))) (* q0y (e/cos (* omega-L t)))))
    :py-expr '(* (* m omega-L) (- (* q0x (e/cos (* omega-L t))) (* q0y (e/sin (* omega-L t)))))}

   ;; Taylor2 with 2nd-order momentum correction:
   ;; q(t) ≈ q₀ + (p₀/m)t − ½·(α·q₀/m)·t²/r₀³
   ;; p(t) ≈ p₀ + F₀·t + ½·(dF/dt)|₀·t²
   ;; where dFx/dt|₀ = α·[(-1/r₀³ + 3q₀x²/r₀⁵)·p₀x + 3q₀x·q₀y·p₀y/r₀⁵] / m
   {:strategy :analytical
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
  ;; Analytical expressions cannot predict full Keplerian trajectories (transcendental).
  ;; We evaluate only the short-horizon portion where Taylor expansion approximations
  ;; (q ≈ q₀ + v₀t + ½·F·t², p ≈ p₀ + F·t) are achievable by GP.
  (try
    (when (analytical-genome-valid? ind)
      (let [{:keys [data]} dataset
            D   (char-scale data)
            fns {:qx (compile-state-fn (:qx-expr ind))
                 :qy (compile-state-fn (:qy-expr ind))
                 :px (compile-state-fn (:px-expr ind))
                 :py (compile-state-fn (:py-expr ind))}
            errors (horizon-errors-analytical fns dataset (short-horizon-times data))]
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
  (or (case (:strategy individual)
        :analytical   (calculate-analytical-fitness   individual dataset)
        :differential (calculate-differential-fitness individual dataset)
        :conserved    (calculate-conserved-fitness    individual dataset)
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

  When :evaluation is :de-driven, scores analytical/conserved genomes by ODE /
  invariant residual on [[phase-states-for-fitness-context]]; differential always 0.

  For :conserved in :data-driven mode, uses a combined per-trajectory + cross-scenario check
  so trivially-constant expressions (e.g. cos(cos(qy/qy)) = const) score 0.

  Optional :aggregate and :percentile (see [[aggregate-scenario-fitness]])."
  [individual datasets & {:keys [aggregate percentile evaluation phase-states]
                          :or {aggregate :min percentile 0.1 evaluation :data-driven}}]
  (if (= (keyword evaluation) :de-driven)
    (or (calculate-de-driven-fitness individual phase-states) 0.0)
    (if (= (:strategy individual) :conserved)
      (or (calculate-conserved-fitness-scenarios individual datasets) 0.0)
      (aggregate-scenario-fitness
       (mapv #(calculate-fitness individual %) datasets)
       :aggregate aggregate
       :percentile percentile))))

(defn- expr-subtrees [expr]
  (let [expr (normalize-expr expr)]
    (if (sequential? expr)
      (cons expr (mapcat expr-subtrees (rest expr)))
      [expr])))

(defn- slotted-subtrees [ind]
  (let [keys (case (:strategy ind)
               :analytical   analytical-expr-keys
               :differential differential-expr-keys
               :conserved    [conserved-expr-key]
               [])]
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

(defn- individual-expr-keys [ind]
  (case (:strategy ind)
    :analytical   analytical-expr-keys
    :differential differential-expr-keys
    :conserved    [conserved-expr-key]
    []))

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
  (+ (* 2.0 (count (filter #(and (sequential? %) (#{'e/sin 'e/cos} (first %)))
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
     :trig? (or (expr-uses-op? expr 'e/sin) (expr-uses-op? expr 'e/cos))
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
  (case slot
    :c-expr (conserved-guess-candidates block-repr)
    (:dqx-expr :dqy-expr :dpx-expr :dpy-expr) (differential-guess-candidates block-repr)
    []))

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

(defn- validate-and-repair
  "Simplify, detect degenerate (constant) sub-expressions, and ensure required
   symbols survive after mutation/crossover."
  [ind]
  (let [strategy (:strategy ind)
        [ks vars] (case strategy
                    :analytical   [analytical-expr-keys   analytical-vars]
                    :differential [differential-expr-keys differential-vars]
                    :conserved    [[conserved-expr-key]   conserved-vars]
                    [nil nil])
        ;; 1. Simplify every expression in the genome; replace any that collapsed
        ;;    to a pure constant (e.g. (- py py) → 0) with a fresh random expr.
        ind' (if ks
               (reduce (fn [acc k]
                         (let [simplified (simplify-expr (get acc k))]
                           (assoc acc k
                                  (if (number? simplified)
                                    (random-expression 3 vars)
                                    simplified))))
                       ind ks)
               ind)
        ;; 2. Re-inject any required symbols that were lost.
    ind'' (sync-block-genome
           (case strategy
             :analytical   (-> ind'
                               (ensure-symbol-coverage ks required-analytical-symbols)
                               ensure-analytical-uses-t)
             :differential (ensure-symbol-coverage ind' ks required-differential-symbols)
             :conserved    (if (conserved-genome-valid? ind')
                             ind'
                             (random-conserved-individual))
             ind'))]
    ind''))

(defn mutate-individual [ind]
  (let [restart-p (if *stagnation-escape?* 0.4 0.2)]
    (-> (cond
          (< (rand) restart-p) (random-individual)
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
  (let [donor-nodes (filterv #(or (coll? %) (symbol? %)) (all-subtrees donor))
        donor-sub   (rand-nth (if (seq donor-nodes) donor-nodes [donor]))]
    (normalize-expr (splice-subtree recipient donor-sub 0.25))))

(defn crossover-individuals
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
                child        (-> {:strategy strategy :expr-blocks child-blocks}
                               (sync-block-genome :from-blocks? true))
                mutated      (reduce (fn [acc k]
                                       (update acc k #(mutate-expr-via-blocks % vars k)))
                                     child ks)]
            (validate-and-repair mutated)))))))


(def default-population-file "data/population.edn")
(def checkpoint-version 8)
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
  (let [ind  (sync-block-genome ind)
        norm (fn [k] [k (normalize-expr (get ind k))])]
    (case (:strategy ind)
      :analytical   (into {:strategy :analytical :expr-blocks (:expr-blocks ind)}
                          (map norm analytical-expr-keys))
      :differential (into {:strategy :differential :expr-blocks (:expr-blocks ind)}
                          (map norm differential-expr-keys))
      :conserved    (into {:strategy :conserved :expr-blocks (:expr-blocks ind)}
                          [(norm conserved-expr-key)])
      (into {:strategy (:strategy ind) :expr-blocks (:expr-blocks ind)}
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
          (when (contains? #{7 8} version)
            {:population (mapv sync-block-genome population)
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
  [{:keys [fresh? seed? path population-size strategy de-driven?]}]
  (let [;; Filter seeds to the active strategy (nil = keep all).
        seeds (when seed?
                (cond->> physics-seeds
                  strategy (filterv #(= (:strategy %) strategy))
                  de-driven? (remove #(= (:strategy %) :differential))))
        base  (if fresh?
                {:population (vec (repeatedly population-size random-individual))
                 :generations-run 0
                 :resumed? false}
                (if-let [{:keys [population generations-run]} (load-population path)]
                  {:population population
                   :generations-run generations-run
                   :resumed? true}
                  {:population (vec (repeatedly population-size random-individual))
                   :generations-run 0
                   :resumed? false}))]
    (if (seq seeds)
      (do (println (str "Injecting " (count seeds) " physics seeds into initial population."))
          (update base :population splice-seeds seeds))
      base)))

(def default-mcts-simulations 64)
(def default-mcts-inject 5)

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
               :mcts? true
               :prompt-each-generation false
               :fitness-mode :random
               :scenario-samples 32
               :fitness-aggregate :min
               :fitness-percentile 10
               :scenario-seed nil
               :guess-mutations? true
               :de-driven? false
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
          "--mcts-until-stop" (recur (assoc opts :mcts-until-stop true) more)
          "--prompt-each-generation" (recur (assoc opts :prompt-each-generation true) more)
          "--fixed-scenarios" (recur (assoc opts :fitness-mode :fixed) more)
          "--random-scenarios" (recur (assoc opts :fitness-mode :random) more)
          "--scenario-samples" (recur (assoc opts :scenario-samples (Long/parseLong (first more))) (rest more))
          "--fitness-aggregate" (recur (assoc opts :fitness-aggregate (keyword (first more))) (rest more))
          "--fitness-percentile" (recur (assoc opts :fitness-percentile (Long/parseLong (first more))) (rest more))
          "--scenario-seed" (recur (assoc opts :scenario-seed (Long/parseLong (first more))) (rest more))
          "--de-driven" (recur (assoc opts :de-driven? true) more)
          "--no-guess" (recur (assoc opts :guess-mutations? false) more)
          "--mcts-simulations" (recur (assoc opts :mcts-simulations (Long/parseLong (first more))) (rest more))
          "--mcts-inject" (recur (assoc opts :mcts-inject (Long/parseLong (first more))) (rest more))
          "--population" (recur (assoc opts :path (first more)) (rest more))
          "--generations" (recur (assoc opts :generations (Long/parseLong (first more))) (rest more))
          "--population-size" (recur (assoc opts :population-size (Long/parseLong (first more))) (rest more))
          "--strategy" (recur (assoc opts :strategy (keyword (first more))) (rest more))
          (throw (ex-info "Unknown argument"
                          {:arg a
                           :hint "--fresh --de-driven --fixed-scenarios --random-scenarios --scenario-samples N --fitness-aggregate min|percentile --fitness-percentile P --scenario-seed N --no-guess --no-mcts ..."})))))))

(defn- distinct-elites
  "Elite tier for one generation. During escape burst, collapse behaviorally identical clones."
  [scored elite-cap probes behavior-diverse? cache]
  (if behavior-diverse?
    (take elite-cap (take-distinct-by-behavior (count scored) scored probes cache))
    (take elite-cap (distinct scored))))

(defn evolve-generation
  [population fitness-ctx population-size generation-index
   & {:keys [extra-immigrants elite-divisor behavior-probes behavior-diverse-elites?
             behavior-cache score-progress? gen-label]
      :or {extra-immigrants 0 elite-divisor 5 behavior-diverse-elites? false
           score-progress? false}}]
  (let [evaluation  (:evaluation fitness-ctx :data-driven)
        de-driven? (= evaluation :de-driven)
        datasets    (when-not de-driven?
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
                     (binding [*fast-differential-fitness?* true
                               *fast-conserved-fitness?* true
                               *fitness-max-states* 64]
                       (map-indexed
                        (fn [i ind]
                          (when score-progress?
                            (do (println (format "    %s scoring %d/%d (%s)..."
                                                 (or gen-label "gen") (inc i) n-pop
                                                 (name (:strategy ind))))
                                (flush)))
                          (assoc ind :fitness
                                 (or (call-with-timeout
                                      (or *fitness-timeout-ms* 0)
                                      (name (:strategy ind))
                                      #(calculate-fitness-scenarios
                                        ind datasets
                                        :aggregate (:aggregate fit-opts)
                                        :percentile (:percentile fit-opts)
                                        :evaluation evaluation
                                        :phase-states phase-states))
                                     0.0)))
                        population)))
        scored      (sort-by :fitness #(compare %2 %1) scored)
        unique-elites (distinct-elites scored elite-cap behavior-probes behavior-diverse-elites? cache)
        n-elites      (count unique-elites)
        breed-slots   (max 0 (- population-size immigrant-n))
        branch        (long (Math/ceil (/ breed-slots (double n-elites))))]
    (vec (take population-size
               (concat
                ;; Elites + their offspring (crossover or mutation).
                ;; Capped at breed-slots so immigrants always get their reserved slots.
                (take breed-slots
                      (mapcat (fn [parent]
                                (cons parent
                                      (take (dec branch)
                                            (repeatedly
                                             #(if (and (> (count unique-elites) 1) (< (rand) 0.7))
                                                ;; Cross-pollination with a randomly chosen
                                                ;; different elite; fall back to mutation if
                                                ;; strategies differ.
                                                (let [other (rand-nth (remove #{parent} unique-elites))]
                                                  (or (crossover-individuals parent other)
                                                      (mutate-individual parent)))
                                                (mutate-individual parent))))))
                              unique-elites))
                ;; Fresh random immigrants — guaranteed to appear every generation.
                (binding [*fast-immigrants?* (or *fast-immigrants?* (> extra-immigrants 5))]
                  (repeatedly immigrant-n random-individual)))))))

(defn- prompt-continue-evolution? []
  (print "  Enter = next generation, q = stop and save: ")
  (flush)
  (not= "q" (clojure.string/trim (or (read-line) ""))))

(defn -main [& args]
  (timbre/merge-config! {:min-level :warn})
  (let [{:keys [fresh? seed? path generations population-size prompt-each-generation
                fitness-context fitness-mode scenario-samples fitness-aggregate fitness-percentile
                strategy guess-mutations?]}
        (parse-args args)
        de-driven? (= :de-driven (:evaluation fitness-context))]
  (binding [*strategy-filter* strategy
            *guess-mutations?* (if (false? guess-mutations?) false *guess-mutations?*)
            *de-driven-search?* de-driven?]
  (when (and de-driven? (= strategy :differential))
    (println "warning: --strategy differential with --de-driven scores 0 (DE is already known)"))
  (let [report-scenarios default-scenarios
        {:keys [population generations-run resumed?]}
        (resolve-initial-population {:fresh? fresh?
                                     :seed?  seed?
                                     :path path
                                     :population-size population-size
                                     :strategy strategy
                                     :de-driven? de-driven?})
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
        ref-datasets (scenarios->datasets default-scenarios)
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
                           pop-with-hof     (if (and @hall-of-fame (not escape?))
                                              (assoc (vec pop) (dec (count pop))
                                                     (:ind @hall-of-fame))
                                              pop)
                           hof-behavior     (when (and escape-burst? @hall-of-fame)
                                              (behavior-key-for (:ind @hall-of-fame)
                                                                behavior-probes-lite
                                                                behavior-cache))
                           pop-diverse      (binding [*fast-immigrants?* (or escape-burst? escape-deep?)]
                                              (cond
                                                escape-deep?
                                                (vec (repeatedly population-size random-individual))

                                                escape-burst?
                                                (let [without-basin (if hof-behavior
                                                                      (vec (remove #(= hof-behavior
                                                                                       (behavior-key-for %
                                                                                                          behavior-probes-lite
                                                                                                          behavior-cache))
                                                                                   pop-with-hof))
                                                                      pop-with-hof)
                                                      n-keep (max 0 (quot population-size 5))
                                                      kept   (vec (take n-keep without-basin))
                                                      n-fresh (min (- population-size (count kept))
                                                                   (quot population-size 3))]
                                                  (into kept (repeatedly n-fresh random-individual)))

                                                :else pop-with-hof))
                           pop' (binding [*stagnation-escape?* (or escape-burst? escape-deep?)
                                          *fast-immigrants?* (or escape-burst? escape-deep?)
                                          *fitness-timeout-ms* (when (or escape-burst? escape-deep?)
                                                                 180000)]
                                  (evolve-generation pop-diverse fitness-context population-size gen-idx
                                                     :extra-immigrants extra-imm
                                                     :elite-divisor elite-divisor
                                                     :behavior-probes behavior-probes-lite
                                                     :behavior-diverse-elites? escape-burst?
                                                     :behavior-cache behavior-cache
                                                     :score-progress? (or escape-burst? escape-deep?)
                                                     :gen-label (format "gen %d" (inc gen-idx))))
                           gens (+ generations-run (inc gen-idx))
                           {:keys [best mean median]} (population-fitness-stats pop')
                           ;; Eval every elite individual on the fixed reference scenarios and
                           ;; pick the best — more robust than trusting training rank alone.
                           elite-inds (filter :fitness pop')
                           eval-pairs (keep (fn [ind]
                                              (let [s (calculate-fitness-scenarios
                                                       ind (when-not de-driven? ref-datasets)
                                                       :aggregate (:aggregate fit-opts)
                                                       :percentile (:percentile fit-opts)
                                                       :evaluation (:evaluation fit-opts)
                                                       :phase-states eval-phase-states)]
                                                (when (pos? s) [ind s])))
                                            elite-inds)
                           [best-eval-ind eval-best] (when (seq eval-pairs)
                                                       (apply max-key second eval-pairs))
                           ;; Update hall of fame if this generation produced a new best.
                           _ (when (and eval-best
                                        (or (nil? @hall-of-fame)
                                            (> eval-best (:eval-fitness @hall-of-fame))))
                               (reset! hall-of-fame {:ind best-eval-ind :eval-fitness eval-best}))
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
        final-datasets (when-not de-driven?
                         (datasets-for-fitness-context fitness-context :generation generations))
        final-probes (build-behavior-probes (or final-datasets ref-datasets))
        ranked (->> final-pop
                    (map (fn [ind]
                           (assoc ind :fitness (calculate-fitness-scenarios
                                                ind (when-not de-driven? final-datasets)
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
                    " | evaluation: DE-driven (analytical ODE residual + conserved ∇C·f)")
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
    (doseq [[i {:keys [strategy fitness] :as ind}] (map-indexed vector top5)]
      (println (str "  #" (inc i) " " (name strategy) " fitness=" (format "%.4f" fitness)))
      (case strategy
        :analytical
        (doseq [k analytical-expr-keys]
          (println (str "    " (name k) ":" (pr-str (normalize-expr (get ind k))))))
        :differential
        (doseq [k differential-expr-keys]
          (println (str "    " (name k) ":" (pr-str (normalize-expr (get ind k))))))
        :conserved
        (println (str "    c-expr:" (pr-str (normalize-expr (get ind conserved-expr-key)))))
        (println "    (unknown strategy)")))
    (println "\nBest on fixed reference scenarios (MSE):")
    (doseq [scenario report-scenarios]
      (let [dataset (scenario-data scenario)
            metrics (evaluate-predictions best dataset)]
        (println (str "  " (name (:id scenario)) " " (name (:strategy metrics))
                      " mse=" (format "%.6f" (:mse metrics))))))))))
