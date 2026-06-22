(ns evophy.coords
  "Phase-space coordinate charts for 2D central-force problems.
  Ground-truth integration stays Cartesian; evolved laws may use the chart
  best suited to each scenario's IC (e.g. polar for axisymmetric orbits).")

(def cartesian-state-vars '[qx qy px py])
(def cartesian-ic-vars '[q0x q0y p0x p0y])
(def cartesian-analytical-expr-keys [:qx-expr :qy-expr :px-expr :py-expr])

(def polar-state-vars '[r theta pr ptheta])
(def polar-ic-vars '[r0 theta0 pr0 ptheta0])
(def polar-analytical-expr-keys [:r-expr :theta-expr :pr-expr :ptheta-expr])

(def charts #{:cartesian :polar})

(defn chart? [x]
  (contains? charts (keyword x)))

(defn preferred-chart
  "Polar when IC has approximate radial symmetry in the x-axis plane (q0y ≈ 0, p0x ≈ 0)."
  [{:keys [q0x q0y p0x p0y]}]
  (let [tol 1e-6]
    (if (and (< (Math/abs (double (or q0y 0.0))) tol)
             (< (Math/abs (double (or p0x 0.0))) tol)
             (pos? (Math/abs (double (or q0x 0.0)))))
      :polar
      :cartesian)))

(defn cartesian->polar
  "Canonical polar state from Cartesian (qx,qy,px,py).
   ptheta is the scalar angular momentum L = qx·py − qy·px."
  [{:keys [qx qy px py]}]
  (let [qx (double qx) qy (double qy) px (double px) py (double py)
        r2 (+ (* qx qx) (* qy qy))
        r  (Math/sqrt (max r2 1e-24))]
    {:r r
     :theta (Math/atan2 qy qx)
     :pr (/ (+ (* qx px) (* qy py)) (max r 1e-24))
     :ptheta (- (* qx py) (* qy px))}))

(defn initial-polar-ics
  "Polar ICs from the first trajectory point (or Cartesian IC map)."
  [data-or-ic]
  (let [{:keys [qx qy px py]} (if (map? data-or-ic)
                                (if (:data data-or-ic)
                                  (first (:data data-or-ic))
                                  data-or-ic)
                                data-or-ic)
        {:keys [r theta pr ptheta]} (cartesian->polar {:qx qx :qy qy :px px :py py})]
    {:r0 r :theta0 theta :pr0 pr :ptheta0 ptheta}))

(defn polar-hamiltonian-deriv
  "Hamiltonian vector field in polar canonical coords (r, θ, p_r, p_θ).
   H = p_r²/(2m) + p_θ²/(2m r²) − α/r."
  [m alpha {:keys [r pr ptheta]}]
  (let [md (double m) ad (double alpha)
        r  (max (double r) 1e-12)
        r2 (* r r) r3 (* r2 r)
        pt (double ptheta)]
    {:dr (/ (double pr) md)
     :dtheta (/ pt (* md r2))
     :dpr (- (/ ad r2) (/ (* pt pt) (* md r3)))
     :dptheta 0.0}))

(defn enrich-state [s]
  (merge s (cartesian->polar s)))

(defn enrich-dataset
  "Attach :chart and polar coords to each trajectory point."
  [dataset]
  (let [chart (preferred-chart dataset)
        data  (mapv enrich-state (:data dataset))]
    (assoc dataset :chart chart :data data)))

(defn law-chart
  "Explicit :chart on the law only. Use [[effective-chart]] during fitness."
  [law]
  (let [c (:chart law)]
    (when (chart? c) (keyword c))))

(defn effective-chart
  "Chart for evaluating a law on a scenario: explicit law :chart, else scenario :chart, else :cartesian."
  [law scenario-chart]
  (or (law-chart law) (when (chart? scenario-chart) (keyword scenario-chart)) :cartesian))

(defn analytical-expr-keys-for-chart [chart]
  (case (keyword chart)
    :polar polar-analytical-expr-keys
    cartesian-analytical-expr-keys))

(defn analytical-exprs-present?
  [law chart]
  (every? #(some? (get law %)) (analytical-expr-keys-for-chart chart)))

(defn dominant-chart
  "Most common :chart among scenario datasets (for immigrants / validation defaults)."
  [datasets]
  (if (empty? datasets)
    :cartesian
    (->> datasets (map #(or (:chart %) :cartesian)) frequencies (apply max-key val))))

(defn state-vars-for-chart [chart]
  (case (keyword chart)
    :polar polar-state-vars
    cartesian-state-vars))

(defn ic-vars-for-chart [chart]
  (case (keyword chart)
    :polar polar-ic-vars
    cartesian-ic-vars))

(defn required-analytical-symbols-for-chart [chart]
  (case (keyword chart)
    :polar '[r0 m]
    '[q0x q0y m]))
