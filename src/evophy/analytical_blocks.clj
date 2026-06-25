(ns evophy.analytical-blocks
  "Catalog of known bound-orbit analytical blocks (circle, ellipse, Taylor, …).
   Emmy validates/simplifies entries at load; injection grafts bound arms into the GP population."
  (:require [evophy.symbolic :as sym]))

(def ^:private cartesian-slots
  [:qx-expr :qy-expr :px-expr :py-expr])

(defn- law [id label tags slots]
  {:id id :label label :regime :bound :tags tags :slots slots})

(def analytical-block-catalog
  "Pre-validated bound analytical laws. Unbound arms are supplied separately (Taylor template)."
  [(law :taylor-bound "First-order Taylor (bound)"
        #{:taylor :polynomial}
        '{:qx-expr (+ q0x (* (e/div p0x m) t))
          :qy-expr (+ q0y (* (e/div p0y m) t))
          :px-expr (+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
          :py-expr (+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))})

   (law :taylor2-bound "Second-order Taylor (bound)"
        #{:taylor :polynomial}
        '{:qx-expr (+ q0x (+ (* (e/div p0x m) t) (* (e/div (* (* -0.5 alpha) q0x) (* m r03)) (* t t))))
          :qy-expr (+ q0y (+ (* (e/div p0y m) t) (* (e/div (* (* -0.5 alpha) q0y) (* m r03)) (* t t))))
          :px-expr (+ p0x (* (e/div (* (* -1.0 alpha) q0x) r03) t))
          :py-expr (+ p0y (* (e/div (* (* -1.0 alpha) q0y) r03) t))})

   (law :circle-omega "Circular orbit R(ωt)·q₀, ω=√(α/mr₀³)"
        #{:circle :trig :closed}
        '{:qx-expr (- (* q0x (e/cos (* omega t))) (* q0y (e/sin (* omega t))))
          :qy-expr (+ (* q0x (e/sin (* omega t))) (* q0y (e/cos (* omega t))))
          :px-expr (* (* -1.0 (* m omega)) (+ (* q0x (e/sin (* omega t))) (* q0y (e/cos (* omega t)))))
          :py-expr (* (* m omega) (- (* q0x (e/cos (* omega t))) (* q0y (e/sin (* omega t)))))})

   (law :ellipse-omega-L "Eccentric orbit R(Ω_L t)·q₀, Ω_L=L/(mr₀²)"
        #{:ellipse :trig :closed}
        '{:qx-expr (- (* q0x (e/cos (* omega-L t))) (* q0y (e/sin (* omega-L t))))
          :qy-expr (+ (* q0x (e/sin (* omega-L t))) (* q0y (e/cos (* omega-L t))))
          :px-expr (* (* -1.0 (* m omega-L)) (+ (* q0x (e/sin (* omega-L t))) (* q0y (e/cos (* omega-L t)))))
          :py-expr (* (* m omega-L) (- (* q0x (e/cos (* omega-L t))) (* q0y (e/sin (* omega-L t)))))})

   (law :ellipse-axis-mix "Ellipse: q along major axis + p/m phase mix"
        #{:ellipse :trig}
        '{:qx-expr (+ (* q0x (e/cos (* omega t))) (* (e/div p0y m) (e/sin (* omega t))))
          :qy-expr (+ (* q0y (e/cos (* omega t))) (* (e/div (- p0x) m) (e/sin (* omega t))))
          :px-expr (+ (* (* -1.0 (* m omega)) q0y (e/sin (* omega t)))
                      (* p0x (e/cos (* omega t))))
          :py-expr (+ (* (* m omega) q0x (e/sin (* omega t)))
                      (* p0y (e/cos (* omega t))))})

   (law :harmonic-shift "Shifted harmonic: q₀ cos ωt + (p₀/mω) sin ωt"
        #{:harmonic :trig}
        '{:qx-expr (+ (* q0x (e/cos (* omega t)))
                      (* (e/div p0x (* m omega)) (e/sin (* omega t))))
          :qy-expr (+ (* q0y (e/cos (* omega t)))
                      (* (e/div p0y (* m omega)) (e/sin (* omega t))))
          :px-expr (+ (* (* -1.0 (* m omega)) q0x (e/sin (* omega t)))
                      (* p0x (e/cos (* omega t))))
          :py-expr (+ (* (* -1.0 (* m omega)) q0y (e/sin (* omega t)))
                      (* p0y (e/cos (* omega t))))})])

(def validated-catalog
  (vec (for [entry analytical-block-catalog
             :when (every? sym/symbolic-validate-expr? (vals (:slots entry)))]
         entry)))

(defn catalog-entry [id]
  (first (filter #(= id (:id %)) validated-catalog)))

(defn catalog-laws
  "Legacy analytical maps for seeds / immigrants."
  []
  (mapv (fn [entry]
          (into {:strategy :analytical :domain :bound}
                (:slots entry)))
        validated-catalog))

(defn random-catalog-entry []
  (rand-nth validated-catalog))

(defn random-catalog-law []
  (into {:strategy :analytical :domain :bound}
        (:slots (random-catalog-entry))))

(defn slot-expr [entry slot]
  (get (:slots entry) slot))

(defn entries-for-slot [slot]
  (filterv #(contains? (:slots %) slot) validated-catalog))

(defn bound-pairs []
  [[:qx-expr :px-expr] [:qy-expr :py-expr]])

(defn pair-slot-updates
  "Map of slot → bound-arm expression for one q/p pair."
  [entry pair]
  (into {}
        (map (fn [slot] [slot (slot-expr entry slot)])
             pair)))

(defn graft-bound-pair
  "Assoc bare bound-arm exprs for one pair (caller wraps e/if + Taylor unbound)."
  [leg entry pair]
  (into leg (pair-slot-updates entry pair)))

(defn graft-catalog-pair [leg & {:keys [entry pair]
                               :or {entry (random-catalog-entry)
                                    pair (rand-nth (bound-pairs))}}]
  (graft-bound-pair leg entry pair))

(defn catalog-block-count []  (count validated-catalog))
