(ns evophy.primitives
  "Graded pure-function extension to the GP expression language.
   Tier 0 = algebra only; higher tiers unlock one Kepler pipeline stage at a time.")

;; ── Complex helpers (unified Kepler oracle only) ─────────────────────────────
;; z = re + i·im.  Hyperbolic sheet: M = −i·M_h, u = i·F, n = (−2E)^(3/2)/(α√m).

(defn- c+ [[ar ai] [br bi]] [(+ ar br) (+ ai bi)])
(defn- c- [[ar ai] [br bi]] [(- ar br) (- ai bi)])
(defn- c* [[ar ai] [br bi]]
  [(- (* ar br) (* ai bi)) (+ (* ar bi) (* ai br))])
(defn- c*sc [s [r i]] [(* s r) (* s i)])
(defn- c-re [[r _]] r)
(defn- c-im [[_ i]] i)

(defn- csin [[r i]]
  [(* (Math/sin r) (Math/cosh i))
   (* (Math/cos r) (Math/sinh i))])

(defn- ccos [[r i]]
  [(* (Math/cos r) (Math/cosh i))
   (- (* (Math/sin r) (Math/sinh i)))])

(defn- cdiv [[ar ai] [br bi]]
  (let [den (+ (* br br) (* bi bi))]
    (if (< den 1e-300)
      [0.0 0.0]
      [(/ (+ (* ar br) (* ai bi)) den)
       (/ (- (* ai br) (* ar bi)) den)])))

(defn- csqrt-principal [[r i]]
  (let [mag (Math/sqrt (+ (* r r) (* i i)))
        ang (Math/atan2 i r)]
    (if (< mag 1e-300)
      [0.0 0.0]
      (let [sm (Math/sqrt mag) half (* 0.5 ang)]
        [(* sm (Math/cos half)) (* sm (Math/sin half))]))))

(defn- mean-motion-complex
  "n = (−2E)^(3/2) / (α√m) — real for E<0, purely imaginary for E>0."
  [E alpha m]
  (let [base (* -2.0 (double E))
        mag  (Math/pow (Math/abs base) 1.5)
        scale (/ 1.0 (* (double alpha) (Math/sqrt (double m))))
        [nr ni] (if (neg? base)
                  [0.0 (- mag)]   ; (−2E)^(3/2) = −i·(2E)^(3/2) for E>0
                  [mag 0.0])]
    (c*sc scale [nr ni])))

(defn- mean-motion-real [E alpha m]
  (/ (Math/pow (* 2.0 (Math/abs (double E))) 1.5)
     (* (double alpha) (Math/sqrt (double m)))))

(defn- M-complex-at-ic [ecc energy u-or-F-at-ic]
  (if (neg? (double energy))
    (let [u (double u-or-F-at-ic)
          M (- u (* (double ecc) (Math/sin u)))]
      [M 0.0])
    (let [F (double u-or-F-at-ic)
          Mh (- (* (double ecc) (Math/sinh F)) F)]
      [0.0 (- Mh)])))

(defn- M-complex-at-t [M0c nc t]
  (c+ M0c (c*sc (double t) nc)))

(defn- wrap-elliptic-M [M]
  (let [twopi (* 2.0 Math/PI)
        M     (double M)]
    (- M (* twopi (Math/floor (/ M twopi))))))

(defn- solve-kepler-unified
  "Invert M = u − e·sin(u) for complex M and u (Newton on ℂ)."
  [ecc Mc]
  (let [e (double ecc)
        Mc (let [[mr mi] Mc]
             (if (< (Math/abs mi) 1e-14)
               [(wrap-elliptic-M mr) 0.0]
               Mc))
        [mr mi] Mc
        u0 (if (< (Math/abs mi) 1e-14)
             [(double mr) 0.0]
             [0.0 (- (double mi))])]
    (loop [u u0 n 0]
      (let [su (csin u)
            cu (ccos u)
            f  (c- (c- u (c*sc e su)) Mc)
            fp (c- [1.0 0.0] (c*sc e cu))
            du (cdiv f fp)]
        (if (or (> n 80)
                (< (+ (Math/abs (c-re du)) (Math/abs (c-im du))) 1e-12))
          u
          (recur (c- u du) (inc n)))))))

(defn- semi-major-a [E alpha]
  (/ (double alpha) (* -2.0 (double E))))

(defn- semi-minor-complex [a ecc]
  (c* [a 0.0] (csqrt-principal [(- 1.0 (* (double ecc) (double ecc))) 0.0])))

(defn- periapsis-from-u [a ecc u n]
  (let [su (csin u)
        cu (ccos u)
        b  (semi-minor-complex a ecc)
        one-ecu (c- [1.0 0.0] (c*sc (double ecc) cu))
        du-dt (cdiv n one-ecu)
        xp (c-re (c* [a 0.0] (c- cu [ecc 0.0])))
        yp (c-re (c* b su))
        vxp (c-re (c* (c*sc (- a) su) du-dt))
        vyp (c-re (c* (c* b cu) du-dt))]
    {:kepler-xp xp :kepler-yp yp :kepler-vxp vxp :kepler-vyp vyp}))

(defn grav2d-energy
  [m alpha qx qy px py]
  (let [r (max (Math/sqrt (+ (* qx qx) (* qy qy))) 1e-12)]
    (+ (/ (+ (* px px) (* py py)) (* 2.0 m))
       (- (/ alpha r)))))

(defn solve-kepler-elliptic
  ^double [^double e ^double M]
  (c-re (solve-kepler-unified e [M 0.0])))

(defn solve-kepler-hyperbolic
  ^double [^double e ^double Mh]
  (c-im (solve-kepler-unified e [0.0 (- Mh)])))

(defn mean-anomaly
  "M(t) = M₀ + n t (real API; hyperbolic uses physical M_h)."
  [n t m0]
  (+ (double m0) (* (double n) (double t))))

(defn anomaly-from-M
  "Invert Kepler's time equation via unified complex oracle."
  [e energy M]
  (let [Mc (if (neg? (double energy))
             [(double M) 0.0]
             [0.0 (- (double M))])
        u  (solve-kepler-unified (double e) Mc)]
    (if (neg? (double energy))
      (c-re u)
      (c-im u))))

(defn- anomaly->u-complex [energy anomaly]
  (if (neg? (double energy))
    [(double anomaly) 0.0]
    [0.0 (double anomaly)]))

(defn periapsis-xp
  [semi-a ecc energy anomaly]
  (let [a-pos (double semi-a)
        signed-a (if (neg? (double energy)) a-pos (- a-pos))
        u (anomaly->u-complex energy anomaly)]
    (c-re (c* [signed-a 0.0] (c- (ccos u) [(double ecc) 0.0])))))

(defn periapsis-yp
  [semi-b bh energy anomaly]
  (let [u (anomaly->u-complex energy anomaly)]
    (if (neg? (double energy))
      (* (double semi-b) (Math/sin (c-re u)))
      (* (double bh) (Math/sinh (c-im u))))))

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
                M0 (c-re (M-complex-at-ic ecc E u0))]
            {:energy E :ecc ecc :cos-om cos-om :sin-om sin-om
             :semi-a a-pos :semi-b b :bh Double/NaN :n-mean n :kepler-M0 M0})
          (let [bh (* a-pos (Math/sqrt (max 0.0 (- (* ecc ecc) 1.0))))
                sh0 (/ y0p bh)
                F0  (Math/log (+ sh0 (Math/sqrt (+ 1.0 (* sh0 sh0)))))
                M0  (- (* ecc (Math/sinh F0)) F0)]
            {:energy E :ecc ecc :cos-om cos-om :sin-om sin-om
             :semi-a a-pos :semi-b Double/NaN :bh bh :n-mean n :kepler-M0 M0}))))))

(defn kepler-orbit-nan []
  {:qx Double/NaN :qy Double/NaN :px Double/NaN :py Double/NaN
   :ecc Double/NaN :cos-om Double/NaN :sin-om Double/NaN
   :semi-a Double/NaN :semi-b Double/NaN :bh Double/NaN :n-mean Double/NaN
   :kepler-u Double/NaN :kepler-F Double/NaN
   :kepler-xp Double/NaN :kepler-yp Double/NaN
   :kepler-vxp Double/NaN :kepler-vyp Double/NaN})

(defn kepler-orbit-at-t
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
      (kepler-orbit-nan)
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
            u        (solve-kepler-unified ecc Mc)
            {:keys [kepler-xp kepler-yp kepler-vxp kepler-vyp]}
            (periapsis-from-u a-signed ecc u nc)
            qx       (- (* kepler-xp cos-om) (* kepler-yp sin-om))
            qy       (+ (* kepler-xp sin-om) (* kepler-yp cos-om))
            vx       (- (* kepler-vxp cos-om) (* kepler-vyp sin-om))
            vy       (+ (* kepler-vxp sin-om) (* kepler-vyp cos-om))
            b-real   (* a-pos (Math/sqrt (max 0.0 (- 1.0 (* ecc ecc)))))]
        (if (neg? E)
          {:qx qx :qy qy :px (* m vx) :py (* m vy)
           :ecc ecc :cos-om cos-om :sin-om sin-om
           :semi-a a-pos :semi-b b-real :bh Double/NaN :n-mean n-real
           :kepler-u (c-re u) :kepler-F Double/NaN
           :kepler-xp kepler-xp :kepler-yp kepler-yp
           :kepler-vxp kepler-vxp :kepler-vyp kepler-vyp}
          (let [bh (* a-pos (Math/sqrt (max 0.0 (- (* ecc ecc) 1.0))))]
            {:qx qx :qy qy :px (* m vx) :py (* m vy)
             :ecc ecc :cos-om cos-om :sin-om sin-om
             :semi-a a-pos :semi-b Double/NaN :bh bh :n-mean n-real
             :kepler-u Double/NaN :kepler-F (c-im u)
             :kepler-xp kepler-xp :kepler-yp kepler-yp
             :kepler-vxp kepler-vxp :kepler-vyp kepler-vyp}))))))

(defn M0-at-ic
  [m alpha q0x q0y p0x p0y]
  (:kepler-M0 (orbital-elements-at-ic m alpha q0x q0y p0x p0y)))

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

(defn expr-primitive-depth
  "Count of graded primitive calls in an expression tree."
  [expr]
  (count (primitive-ops-in-expr expr)))

(defn expr-primitive-valid?
  [expr unlocked-tier]
  (let [ops (primitive-ops-in-expr expr)]
    (and (<= (count ops) (max-pipeline-depth-for-tier unlocked-tier))
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
  [op]
  (case op
    e/mean-anomaly '[n-mean t kepler-M0]
    e/anomaly      '[ecc energy kepler-M]
    e/periapsis-xp '[semi-a ecc energy kepler-u]
    e/periapsis-yp '[semi-b bh energy kepler-u]
    e/lab-qx       '[cos-om sin-om kepler-xp kepler-yp]
    e/lab-qy       '[sin-om cos-om kepler-xp kepler-yp]
    []))

(defn random-primitive-expr
  "One graded primitive call; optional single nested lower-tier arg (pipeline growth)."
  [unlocked-tier vars depth]
  (when (pos? unlocked-tier)
    (let [ops   (unlocked-primitive-ops unlocked-tier)
          op    (rand-nth ops)
          arity (primitive-arity op)
          fill  (fn []
                  (if (and (pos? depth)
                           (< (rand) 0.35)
                           (pos? unlocked-tier))
                    (or (random-primitive-expr (dec unlocked-tier) vars (dec depth))
                        (rand-nth vars))
                    (rand-nth (concat vars (default-arg-symbols op)))))
          args  (vec (repeatedly arity fill))]
      (cons op args))))

(defn kepler-pipeline-expr
  "Tier-5 reference: composed M → anomaly → periapsis → lab (catalog / seed)."
  [slot]
  (let [M* (list 'e/mean-anomaly 'n-mean 't 'kepler-M0)
        u* (list 'e/anomaly 'ecc 'energy M*)]
    (case slot
      :qx-expr (list 'e/lab-qx 'cos-om 'sin-om
                     (list 'e/periapsis-xp 'semi-a 'ecc 'energy u*)
                     (list 'e/periapsis-yp 'semi-b 'bh 'energy 'kepler-u))
      :qy-expr (list 'e/lab-qy 'sin-om 'cos-om
                     (list 'e/periapsis-xp 'semi-a 'ecc 'energy u*)
                     (list 'e/periapsis-yp 'semi-b 'bh 'energy 'kepler-u))
      :px-expr '(* m (- (* kepler-vxp cos-om) (* kepler-vyp sin-om)))
      :py-expr '(* m (+ (* kepler-vxp sin-om) (* kepler-vyp cos-om)))
      M*)))
