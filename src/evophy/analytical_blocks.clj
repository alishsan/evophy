(ns evophy.analytical-blocks
  "Catalog of known analytical blocks (bound orbits, hyperbolic flybys, …).
   Emmy validates/simplifies entries at load; injection grafts arms into the GP population."
  (:require [evophy.symbolic :as sym]))

(def ^:private cartesian-slots
  [:qx-expr :qy-expr :px-expr :py-expr])

(defn- law [id label regime tags slots]
  {:id id :label label :regime regime :tags tags :slots slots})

(def taylor-bound-slots
  "First-order linear bound arms — internal default when wrapping unbound catalog entries (not in catalog)."
  '{:qx-expr (+ q0x (* (e/div p0x m) t))
    :qy-expr (+ q0y (* (e/div p0y m) t))
    :px-expr (+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
    :py-expr (+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))})

(def analytical-block-catalog
  "Pre-validated analytical laws tagged :bound or :unbound (no Taylor expansions)."
  [(law :circle-omega "Circular orbit R(ωt)·q₀, ω=√(α/mr₀³)"
        :bound #{:circle :trig :closed}
        '{:qx-expr (- (* q0x (e/cos (* omega t))) (* q0y (e/sin (* omega t))))
          :qy-expr (+ (* q0x (e/sin (* omega t))) (* q0y (e/cos (* omega t))))
          :px-expr (* (* -1.0 (* m omega)) (+ (* q0x (e/sin (* omega t))) (* q0y (e/cos (* omega t)))))
          :py-expr (* (* m omega) (- (* q0x (e/cos (* omega t))) (* q0y (e/sin (* omega t)))))})

   (law :ellipse-omega-L "Eccentric orbit R(Ω_L t)·q₀, Ω_L=L/(mr₀²)"
        :bound #{:ellipse :trig :closed}
        '{:qx-expr (- (* q0x (e/cos (* omega-L t))) (* q0y (e/sin (* omega-L t))))
          :qy-expr (+ (* q0x (e/sin (* omega-L t))) (* q0y (e/cos (* omega-L t))))
          :px-expr (* (* -1.0 (* m omega-L)) (+ (* q0x (e/sin (* omega-L t))) (* q0y (e/cos (* omega-L t)))))
          :py-expr (* (* m omega-L) (- (* q0x (e/cos (* omega-L t))) (* q0y (e/sin (* omega-L t)))))})

   (law :ellipse-axis-mix "Ellipse: q along major axis + p/m phase mix"
        :bound #{:ellipse :trig}
        '{:qx-expr (+ (* q0x (e/cos (* omega t))) (* (e/div p0y m) (e/sin (* omega t))))
          :qy-expr (+ (* q0y (e/cos (* omega t))) (* (e/div (- p0x) m) (e/sin (* omega t))))
          :px-expr (+ (* (* -1.0 (* m omega)) q0y (e/sin (* omega t)))
                      (* p0x (e/cos (* omega t))))
          :py-expr (+ (* (* m omega) q0x (e/sin (* omega t)))
                      (* p0y (e/cos (* omega t))))})

   (law :harmonic-shift "Shifted harmonic: q₀ cos ωt + (p₀/mω) sin ωt"
        :bound #{:harmonic :trig}
        '{:qx-expr (+ (* q0x (e/cos (* omega t)))
                      (* (e/div p0x (* m omega)) (e/sin (* omega t))))
          :qy-expr (+ (* q0y (e/cos (* omega t)))
                      (* (e/div p0y (* m omega)) (e/sin (* omega t))))
          :px-expr (+ (* (* -1.0 (* m omega)) q0x (e/sin (* omega t)))
                      (* p0x (e/cos (* omega t))))
          :py-expr (+ (* (* -1.0 (* m omega)) q0y (e/sin (* omega t)))
                      (* p0y (e/cos (* omega t))))})

   (law :ellipse-conic "Elliptic conic: a(cos u−e), b sin u rotated by ω"
        :bound #{:ellipse :conic :trig :closed}
        '{:qx-expr (- (* cos-om (* semi-a (- (e/cos kepler-u) ecc)))
                      (* sin-om (* semi-b (e/sin kepler-u))))
          :qy-expr (+ (* sin-om (* semi-a (- (e/cos kepler-u) ecc)))
                      (* cos-om (* semi-b (e/sin kepler-u))))
          :px-expr (* m (- (* kepler-vxp cos-om) (* kepler-vyp sin-om)))
          :py-expr (* m (+ (* kepler-vxp sin-om) (* kepler-vyp cos-om)))})

   (law :hyperbola-conic "Hyperbolic conic: a(e−cosh F), b_h sinh F rotated by ω"
        :unbound #{:hyperbola :conic :sinh :open}
        '{:qx-expr (- (* cos-om (* semi-a (- ecc (e/cosh kepler-F))))
                      (* sin-om (* bh (e/sinh kepler-F))))
          :qy-expr (+ (* sin-om (* semi-a (- ecc (e/cosh kepler-F))))
                      (* cos-om (* bh (e/sinh kepler-F))))
          :px-expr (* m (- (* kepler-vxp cos-om) (* kepler-vyp sin-om)))
          :py-expr (* m (+ (* kepler-vxp sin-om) (* kepler-vyp cos-om)))})])

(def ^:private taylor-unbound-slots
  "First-order Taylor unbound arms — legacy template default."
  '{:qx-expr (+ q0x (* (e/div p0x m) t))
    :qy-expr (+ q0y (* (e/div p0y m) t))
    :px-expr (+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
    :py-expr (+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))})

(def validated-catalog
  (vec (for [entry analytical-block-catalog
             :when (every? sym/symbolic-validate-expr? (vals (:slots entry)))]
         entry)))

(defn catalog-entry [id]
  (first (filter #(= id (:id %)) validated-catalog)))

(defn hyperbola-conic-slots
  "Unbound-arm slot map for :hyperbola-conic (nil when entry missing)."
  []
  (:slots (catalog-entry :hyperbola-conic)))

(defn kepler-conic-entry
  "Synthetic both-regimes entry: :ellipse-conic bound + :hyperbola-conic unbound."
  []
  (when-let [ell (catalog-entry :ellipse-conic)]
    (when-let [hyp (catalog-entry :hyperbola-conic)]
      {:id :kepler-conic
       :label "Kepler conic: elliptic (sin/cos u) bound + hyperbolic (sinh/cosh F) unbound"
       :regime :any
       :tags #{:kepler :conic}
       :bound-slots (:slots ell)
       :unbound-slots (:slots hyp)})))

(defn kepler-conic-valid?
  []
  (boolean (kepler-conic-entry)))

(defn regime->domain [regime]
  (case regime
    :unbound :unbound
    :bound :bound
    :any))

(defn catalog-laws
  "Legacy analytical maps for seeds / immigrants."
  []
  (mapv (fn [entry]
          (into {:strategy :analytical :domain (regime->domain (:regime entry))}
                (:slots entry)))
        validated-catalog))

(defn random-catalog-entry []
  (rand-nth validated-catalog))

(defn random-catalog-law []
  (let [entry (random-catalog-entry)]
    (into {:strategy :analytical :domain (regime->domain (:regime entry))}
          (:slots entry))))

(defn slot-expr [entry slot]
  (get (:slots entry) slot))

(defn entries-for-slot [slot]
  (filterv #(contains? (:slots %) slot) validated-catalog))

(defn entries-for-regime [regime]
  (filterv #(= regime (:regime %)) validated-catalog))

(defn bound-pairs []
  [[:qx-expr :px-expr] [:qy-expr :py-expr]])

(defn unbound-pairs []
  (bound-pairs))

(defn pair-slot-updates
  "Map of slot → arm expression for one q/p pair."
  [entry pair]
  (into {}
        (map (fn [slot] [slot (slot-expr entry slot)])
             pair)))

(defn graft-bound-pair
  "Assoc bare bound-arm exprs for one pair (caller wraps e/if + linear unbound template)."
  [leg entry pair]
  (into leg (pair-slot-updates entry pair)))

(defn graft-unbound-pair
  "Assoc bare unbound-arm exprs for one pair (caller wraps e/if + linear bound template)."
  [leg entry pair]
  (graft-bound-pair leg entry pair))

(defn graft-catalog-pair [leg & {:keys [entry pair]
                               :or {entry (random-catalog-entry)
                                    pair (rand-nth (bound-pairs))}}]
  (if (= :unbound (:regime entry))
    (graft-unbound-pair leg entry pair)
    (graft-bound-pair leg entry pair)))

(defn catalog-block-count []  (count validated-catalog))

(defn injection-catalog-entries
  "Validated single-regime entries plus optional synthetic :kepler-conic."
  []
  (if-let [kc (kepler-conic-entry)]
    (conj validated-catalog kc)
    validated-catalog))
