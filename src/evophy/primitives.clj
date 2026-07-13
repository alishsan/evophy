(ns evophy.primitives
  "Graded pure-function extension to the GP expression language.
   Tier 0 = algebra only; higher tiers unlock one orbital pipeline stage at a time."
  (:require [emmy.complex :as c]
            [emmy.generic :as g]))

;; ── Unified conic oracle (Emmy complex + Java Math) ─────────────────────────

(defn- ->cx ^org.apache.commons.math3.complex.Complex [re im]
  (c/complex (double re) (double im)))

(defn- wrap-elliptic-M [M]
  (let [twopi (* 2.0 Math/PI)
        M     (double M)]
    (- M (* twopi (Math/floor (/ M twopi))))))

(defn- mean-motion-complex
  "n = (−2E)^(3/2) / (α√m) — real for E<0, purely imaginary for E>0."
  [E alpha m]
  (let [scale (/ 1.0 (* (double alpha) (Math/sqrt (double m))))]
    (g/* (c/complex scale 0.0)
         (g/expt (c/complex (* -2.0 (double E)) 0.0) 1.5))))

(defn- mean-motion-real [E alpha m]
  (/ (Math/pow (* 2.0 (Math/abs (double E))) 1.5)
     (* (double alpha) (Math/sqrt (double m)))))

(defn- M-complex-at-ic [ecc energy u-or-F-at-ic]
  (if (neg? (double energy))
    (let [u (double u-or-F-at-ic)]
      (->cx (- u (* (double ecc) (Math/sin u))) 0.0))
    (let [F (double u-or-F-at-ic)
          Mh (- (* (double ecc) (Math/sinh F)) F)]
      (->cx 0.0 (- Mh)))))

(defn- M-complex-at-t [M0c nc t]
  (g/+ M0c (g/* (c/complex (double t) 0.0) nc)))

(defn- solve-anomaly-unified
  "Invert M = u − e·sin(u) for complex M and u (Newton; Emmy g/sin, g/cos)."
  [ecc Mc]
  (let [e (double ecc)
        Mc (if (< (Math/abs (c/imaginary Mc)) 1e-14)
             (->cx (wrap-elliptic-M (c/real Mc)) 0.0)
             Mc)
        u0 (if (< (Math/abs (c/imaginary Mc)) 1e-14)
             (->cx (c/real Mc) 0.0)
             (->cx 0.0 (- (c/imaginary Mc))))]
    (loop [u u0 n 0]
      (let [su (g/sin u)
            cu (g/cos u)
            f  (g/- (g/- u (g/* e su)) Mc)
            fp (g/- (c/complex 1 0) (g/* e cu))
            du (g// f fp)]
        (if (or (> n 80)
                (< (+ (Math/abs (c/real du)) (Math/abs (c/imaginary du))) 1e-12))
          u
          (recur (g/- u du) (inc n)))))))

(defn- semi-major-a [E alpha]
  (/ (double alpha) (* -2.0 (double E))))

(defn- semi-minor-complex [a ecc]
  (g/* (c/complex (double a) 0.0)
       (g/sqrt (c/complex (- 1.0 (* (double ecc) (double ecc))) 0.0))))

(defn- periapsis-from-u [a ecc u nc]
  (let [su (g/sin u)
        cu (g/cos u)
        b  (semi-minor-complex a ecc)
        one-ecu (g/- (c/complex 1 0) (g/* (double ecc) cu))
        du-dt (g// nc one-ecu)
        xp (c/real (g/* (c/complex (double a) 0.0)
                        (g/- cu (c/complex (double ecc) 0.0))))
        yp (c/real (g/* b su))
        vxp (c/real (g/* (g/* (c/complex (- (double a)) 0.0) su) du-dt))
        vyp (c/real (g/* (g/* b cu) du-dt))]
    {:peri-xp xp :peri-yp yp :peri-vxp vxp :peri-vyp vyp}))

(defn grav2d-energy
  [m alpha qx qy px py]
  (let [r (max (Math/sqrt (+ (* qx qx) (* qy qy))) 1e-12)]
    (+ (/ (+ (* px px) (* py py)) (* 2.0 m))
       (- (/ alpha r)))))

(defn mean-anomaly
  "M(t) = M₀ + n t (real API; hyperbolic uses physical M_h)."
  [n t m0]
  (+ (double m0) (* (double n) (double t))))

(defn anomaly-from-M
  "Invert the time equation via unified complex oracle (Emmy)."
  [e energy M]
  (let [Mc (if (neg? (double energy))
             (->cx (double M) 0.0)
             (->cx 0.0 (- (double M))))
        u  (solve-anomaly-unified (double e) Mc)]
    (if (neg? (double energy))
      (c/real u)
      (c/imaginary u))))

(defn- anomaly->cx [energy anomaly]
  (if (neg? (double energy))
    (->cx (double anomaly) 0.0)
    (->cx 0.0 (double anomaly))))

(defn periapsis-xp
  [semi-a ecc energy anomaly]
  (let [a-pos (double semi-a)
        a-signed (if (neg? (double energy)) a-pos (- a-pos))
        u (anomaly->cx energy anomaly)]
    (c/real (g/* (c/complex a-signed 0.0)
                 (g/- (g/cos u) (c/complex (double ecc) 0.0))))))

(defn periapsis-yp
  [semi-b bh energy anomaly]
  (let [u (double anomaly)]
    (if (neg? (double energy))
      (* (double semi-b) (Math/sin u))
      (* (double bh) (Math/sinh u)))))

(defn lab-qx
  [cos-om sin-om xp yp]
  (- (* (double cos-om) (double xp))
     (* (double sin-om) (double yp))))

(defn lab-qy
  [sin-om cos-om xp yp]
  (+ (* (double sin-om) (double xp))
     (* (double cos-om) (double yp))))

(defn orbital-elements-at-ic
  "Orbital elements and mean anomaly M₀ at t=0 from initial conditions."
  [m alpha q0x q0y p0x p0y]
  (let [q0x (double q0x) q0y (double q0y)
        p0x (double p0x) p0y (double p0y)
        m   (double m)   alpha (double alpha)
        r0  (Math/sqrt (+ (* q0x q0x) (* q0y q0y)))
        L   (- (* q0x p0y) (* q0y p0x))
        E   (grav2d-energy m alpha q0x q0y p0x p0y)
        ex  (- (/ (* p0y L) (* m alpha)) (/ q0x r0))
        ey  (- (/ (* (- p0x) L) (* m alpha)) (/ q0y r0))
        ecc (Math/sqrt (+ (* ex ex) (* ey ey)))]
    (if (or (<= (Math/abs L) 1e-8)
            (<= (Math/abs (- ecc 1.0)) 1e-5))
      nil
      (let [omega  (Math/atan2 ey ex)
            cos-om (Math/cos omega)
            sin-om (Math/sin omega)
            x0p    (+ (* q0x cos-om) (* q0y sin-om))
            y0p    (+ (* (- q0x) sin-om) (* q0y cos-om))
            a-pos  (Math/abs (semi-major-a E alpha))
            n      (mean-motion-real E alpha m)]
        (if (neg? E)
          (let [b  (* a-pos (Math/sqrt (max 0.0 (- 1.0 (* ecc ecc)))))
                u0 (Math/atan2 (/ y0p b) (+ (/ x0p a-pos) ecc))
                M0 (c/real (M-complex-at-ic ecc E u0))]
            {:energy E :ecc ecc :cos-om cos-om :sin-om sin-om
             :semi-a a-pos :semi-b b :bh Double/NaN :n-mean n :mean-M0 M0})
          (let [bh (* a-pos (Math/sqrt (max 0.0 (- (* ecc ecc) 1.0))))
                sh0 (/ y0p bh)
                F0  (Math/log (+ sh0 (Math/sqrt (+ 1.0 (* sh0 sh0)))))
                M0  (- (* ecc (Math/sinh F0)) F0)]
            {:energy E :ecc ecc :cos-om cos-om :sin-om sin-om
             :semi-a a-pos :semi-b Double/NaN :bh bh :n-mean n :mean-M0 M0}))))))

(defn orbit-nan []
  {:qx Double/NaN :qy Double/NaN :px Double/NaN :py Double/NaN
   :ecc Double/NaN :cos-om Double/NaN :sin-om Double/NaN
   :semi-a Double/NaN :semi-b Double/NaN :bh Double/NaN :n-mean Double/NaN
   :ecc-anom Double/NaN :hyp-anom Double/NaN
   :peri-xp Double/NaN :peri-yp Double/NaN
   :peri-vxp Double/NaN :peri-vyp Double/NaN})

(defn orbit-at-t
  "Exact Kepler conic via unified complex oracle: M = u − e·sin u, n = (−2E)^(3/2)/(α√m)."
  [m alpha q0x q0y p0x p0y t]
  (let [q0x (double q0x) q0y (double q0y)
        p0x (double p0x) p0y (double p0y)
        m   (double m)   alpha (double alpha) t (double t)
        r0  (Math/sqrt (+ (* q0x q0x) (* q0y q0y)))
        L   (- (* q0x p0y) (* q0y p0x))
        E   (grav2d-energy m alpha q0x q0y p0x p0y)
        ex  (- (/ (* p0y L) (* m alpha)) (/ q0x r0))
        ey  (- (/ (* (- p0x) L) (* m alpha)) (/ q0y r0))
        ecc (Math/sqrt (+ (* ex ex) (* ey ey)))]
    (if (or (<= (Math/abs L) 1e-8)
            (<= (Math/abs (- ecc 1.0)) 1e-5))
      (orbit-nan)
      (let [omega  (Math/atan2 ey ex)
            cos-om (Math/cos omega)
            sin-om (Math/sin omega)
            x0p    (+ (* q0x cos-om) (* q0y sin-om))
            y0p    (+ (* (- q0x) sin-om) (* q0y cos-om))
            a-signed (semi-major-a E alpha)
            a-pos    (Math/abs a-signed)
            nc       (mean-motion-complex E alpha m)
            n-real   (mean-motion-real E alpha m)
            u0       (if (neg? E)
                       (let [b (* a-pos (Math/sqrt (max 0.0 (- 1.0 (* ecc ecc)))))]
                         (Math/atan2 (/ y0p b) (+ (/ x0p a-pos) ecc)))
                       (let [bh (* a-pos (Math/sqrt (max 0.0 (- (* ecc ecc) 1.0))))
                             sh0 (/ y0p bh)]
                         (Math/log (+ sh0 (Math/sqrt (+ 1.0 (* sh0 sh0)))))))
            M0c      (M-complex-at-ic ecc E u0)
            Mc       (M-complex-at-t M0c nc t)
            u        (solve-anomaly-unified ecc Mc)
            {:keys [peri-xp peri-yp peri-vxp peri-vyp]}
            (periapsis-from-u a-signed ecc u nc)
            qx       (- (* peri-xp cos-om) (* peri-yp sin-om))
            qy       (+ (* peri-xp sin-om) (* peri-yp cos-om))
            vx       (- (* peri-vxp cos-om) (* peri-vyp sin-om))
            vy       (+ (* peri-vxp sin-om) (* peri-vyp cos-om))
            b-real   (* a-pos (Math/sqrt (max 0.0 (- 1.0 (* ecc ecc)))))]
        (if (neg? E)
          {:qx qx :qy qy :px (* m vx) :py (* m vy)
           :ecc ecc :cos-om cos-om :sin-om sin-om
           :semi-a a-pos :semi-b b-real :bh Double/NaN :n-mean n-real
           :ecc-anom (c/real u) :hyp-anom Double/NaN
           :peri-xp peri-xp :peri-yp peri-yp
           :peri-vxp peri-vxp :peri-vyp peri-vyp}
          (let [bh (* a-pos (Math/sqrt (max 0.0 (- (* ecc ecc) 1.0))))]
            {:qx qx :qy qy :px (* m vx) :py (* m vy)
             :ecc ecc :cos-om cos-om :sin-om sin-om
             :semi-a a-pos :semi-b Double/NaN :bh bh :n-mean n-real
             :ecc-anom Double/NaN :hyp-anom (c/imaginary u)
             :peri-xp peri-xp :peri-yp peri-yp
             :peri-vxp peri-vxp :peri-vyp peri-vyp}))))))

(defn M0-at-ic
  [m alpha q0x q0y p0x p0y]
  (:mean-M0 (orbital-elements-at-ic m alpha q0x q0y p0x p0y)))

;; ── Graded GP registry ────────────────────────────────────────────────────────

(def primitive-specs
  "Tier N unlocks all primitives with :tier ≤ N. Pipeline depth capped separately."
  '{e/mean-anomaly {:tier 1 :arity 3 :impl mean-anomaly}
    e/anomaly      {:tier 2 :arity 3 :impl anomaly-from-M}
    e/periapsis-xp {:tier 3 :arity 4 :impl periapsis-xp}
    e/periapsis-yp {:tier 4 :arity 4 :impl periapsis-yp}
    e/lab-qx       {:tier 5 :arity 4 :impl lab-qx}
    e/lab-qy       {:tier 5 :arity 4 :impl lab-qy}})

(def max-primitive-pipeline-depth 5)

(def oracle-anomaly-leaf-syms
  "GP-forbidden in pipeline-only mode — anomalies must come from e/anomaly."
  '#{ecc-anom hyp-anom})

(def oracle-position-leaf-syms
  "GP-forbidden in pipeline-only q slots — position must come from e/periapsis-* / e/lab-*."
  '#{peri-xp peri-yp})

(def oracle-mean-leaf-syms
  "GP-forbidden in pipeline-only mode — use e/mean-anomaly instead of mean-M."
  '#{mean-M})

(defn- expr-mentions-sym? [expr sym]
  (cond
    (= expr sym) true
    (and (sequential? expr) (seq expr))
    (some #(expr-mentions-sym? % sym) (rest expr))
    :else false))

(defn expr-uses-forbidden-oracle-leaf?
  [expr forbidden-syms]
  (boolean (some #(expr-mentions-sym? expr %) forbidden-syms)))

(defn pipeline-only-forbidden-syms
  [slot-key]
  (into oracle-anomaly-leaf-syms
        (cond-> oracle-mean-leaf-syms
          (#{:qx-expr :qy-expr} slot-key) (into oracle-position-leaf-syms))))

(defn primitive-op? [op]
  (contains? primitive-specs op))

(defn primitive-tier [op]
  (get-in primitive-specs [op :tier] 0))

(defn primitive-arity [op]
  (get-in primitive-specs [op :arity]))

(defn unlocked-primitive-ops
  [tier]
  (vec (for [[op spec] primitive-specs
             :when (<= (:tier spec) tier)]
         op)))

(defn max-pipeline-depth-for-tier
  [tier]
  (long (max 0 (min tier max-primitive-pipeline-depth))))

(defn primitive-tier-for-hof
  "Auto-unlock tiers from hall-of-fame fitness (gradual, one stage at a time)."
  [hof]
  (cond
    (>= hof 0.98) 5
    (>= hof 0.90) 4
    (>= hof 0.75) 3
    (>= hof 0.50) 2
    (>= hof 0.25) 1
    :else 0))

(defn- primitive-ops-in-expr
  [expr]
  (filter primitive-op?
          (map first (filter sequential? (tree-seq sequential? seq expr)))))

(defn- expr-primitive-path-depths
  "Max primitive ops on any root-to-leaf path (shared subtrees counted per path)."
  [expr]
  (cond
    (not (sequential? expr)) [0]
    (primitive-op? (first expr))
    (let [child-depths (when (seq (rest expr))
                         (mapcat expr-primitive-path-depths (rest expr)))
          base (if (seq child-depths) (apply max child-depths) 0)]
      [(inc base)])
    :else (if (seq (rest expr))
            (mapcat expr-primitive-path-depths (rest expr))
            [0])))

(defn expr-primitive-depth
  "Max graded-primitive pipeline depth on any root-to-leaf path."
  [expr]
  (apply max 0 (expr-primitive-path-depths expr)))

(defn expr-primitive-valid?
  [expr unlocked-tier]
  (let [ops (primitive-ops-in-expr expr)]
    (and (<= (expr-primitive-depth expr) (max-pipeline-depth-for-tier unlocked-tier))
         (every? #(<= (primitive-tier %) unlocked-tier) ops))))

(defn rewrite-primitives-in-expr
  "Expand e/primitive ops to evophy.core eval names (mean-anomaly, …)."
  [expr]
  (clojure.walk/postwalk
   (fn [x]
     (if (and (sequential? x) (seq x))
       (let [[op & args] x]
         (if-some [impl (get-in primitive-specs [op :impl])]
           (cons impl args)
           x))
       x))
   expr))

(defn- default-arg-symbols
  [op & {:keys [pipeline-only?]}]
  (case op
    e/mean-anomaly '[n-mean t mean-M0]
    e/anomaly      (if pipeline-only?
                     '[ecc energy n-mean]
                     '[ecc energy mean-M])
    e/periapsis-xp (if pipeline-only?
                     '[semi-a ecc energy n-mean]
                     '[semi-a ecc energy ecc-anom])
    e/periapsis-yp (if pipeline-only?
                     '[semi-b bh energy n-mean]
                     '[semi-b bh energy ecc-anom])
    e/lab-qx       (if pipeline-only?
                     '[cos-om sin-om semi-a ecc]
                     '[cos-om sin-om peri-xp peri-yp])
    e/lab-qy       (if pipeline-only?
                     '[sin-om cos-om semi-b bh]
                     '[sin-om cos-om peri-xp peri-yp])
    []))

(defn random-primitive-expr
  "One graded primitive call; optional single nested lower-tier arg (pipeline growth)."
  [unlocked-tier vars depth & {:keys [pipeline-only?]}]
  (when (pos? unlocked-tier)
    (let [ops   (unlocked-primitive-ops unlocked-tier)
          op    (rand-nth ops)
          arity (primitive-arity op)
          fill  (fn []
                  (if (and (pos? depth)
                           (< (rand) 0.35)
                           (pos? unlocked-tier))
                    (or (random-primitive-expr (dec unlocked-tier) vars (dec depth)
                                               :pipeline-only? pipeline-only?)
                        (rand-nth vars))
                    (rand-nth (concat vars (default-arg-symbols op :pipeline-only? pipeline-only?)))))
          args  (vec (repeatedly arity fill))]
      (cons op args))))

(defn graded-pipeline-expr
  "Tier-5 reference: composed M → anomaly → periapsis → lab (catalog / seed)."
  [slot]
  (let [M* (list 'e/mean-anomaly 'n-mean 't 'mean-M0)
        u* (list 'e/anomaly 'ecc 'energy M*)]
    (case slot
      :qx-expr (list 'e/lab-qx 'cos-om 'sin-om
                     (list 'e/periapsis-xp 'semi-a 'ecc 'energy u*)
                     (list 'e/periapsis-yp 'semi-b 'bh 'energy u*))
      :qy-expr (list 'e/lab-qy 'sin-om 'cos-om
                     (list 'e/periapsis-xp 'semi-a 'ecc 'energy u*)
                     (list 'e/periapsis-yp 'semi-b 'bh 'energy u*))
      :px-expr '(* m (- (* peri-vxp cos-om) (* peri-vyp sin-om)))
      :py-expr '(* m (+ (* peri-vxp sin-om) (* peri-vyp cos-om)))
      M*)))
