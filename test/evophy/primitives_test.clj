(ns evophy.primitives-test
  (:require [clojure.test :refer :all]
            [evophy.core :as core]
            [evophy.primitives :as prims]))

(deftest mean-anomaly-and-anomaly-match-ecc-anom
  (let [m 1.0 a 1.0 q0x 1.2 q0y 0.0 p0x 0.0 p0y 0.9 t 0.5
        el (prims/orbital-elements-at-ic m a q0x q0y p0x p0y)
        M  (prims/mean-anomaly (:n-mean el) t (:mean-M0 el))
        u  (prims/anomaly-from-M (:ecc el) (:energy el) M)
        ko (#'core/orbit-at-t m a q0x q0y p0x p0y t)]
    (is (some? el))
    (is (some? ko))
    (is (< (Math/abs (- u (:ecc-anom ko))) 1e-9))))

(deftest periapsis-pipeline-matches-orbit-position
  (let [m 1.0 a 1.0 q0x 1.2 q0y 0.0 p0x 0.0 p0y 0.9 t 0.5
        el (prims/orbital-elements-at-ic m a q0x q0y p0x p0y)
        M  (prims/mean-anomaly (:n-mean el) t (:mean-M0 el))
        u  (prims/anomaly-from-M (:ecc el) (:energy el) M)
        xp (prims/periapsis-xp (:semi-a el) (:ecc el) (:energy el) u)
        yp (prims/periapsis-yp (:semi-b el) (:bh el) (:energy el) u)
        ko (#'core/orbit-at-t m a q0x q0y p0x p0y t)]
    (is (< (Math/abs (- xp (:peri-xp ko))) 1e-9))
    (is (< (Math/abs (- yp (:peri-yp ko))) 1e-9))))

(deftest graded-pipeline-depth-capped
  (let [deep '(e/lab-qx cos-om sin-om
              (e/periapsis-xp semi-a ecc energy
                (e/anomaly ecc energy
                  (e/mean-anomaly n-mean t mean-M0)))
              (e/periapsis-yp semi-b bh energy
                (e/anomaly ecc energy
                  (e/mean-anomaly n-mean t mean-M0))))]
    (is (prims/expr-primitive-valid? deep 5))
    (is (not (prims/expr-primitive-valid? deep 2)))))

(deftest pipeline-only-rejects-oracle-leaves
  (let [cheat-qx '(- (* cos-om (* semi-a (- (e/cos ecc-anom) ecc)))
                      (* sin-om (* semi-b (e/sin ecc-anom))))
        pipeline (prims/graded-pipeline-expr :qx-expr)]
    (is (prims/expr-uses-forbidden-oracle-leaf?
         cheat-qx (prims/pipeline-only-forbidden-syms :qx-expr)))
    (is (not (prims/expr-uses-forbidden-oracle-leaf?
              pipeline (prims/pipeline-only-forbidden-syms :qx-expr))))))

(deftest pipeline-only-genome-valid-for-graded-pipeline
  (binding [core/*primitive-tier* 5
            core/*primitive-pipeline-only?* true
            core/*both-regimes?* false]
    (let [expr (prims/graded-pipeline-expr :qx-expr)
          ind {:strategy :analytical :domain :bound
               :qx-expr expr :qy-expr expr
               :px-expr '(* m (- (* peri-vxp cos-om) (* peri-vyp sin-om)))
               :py-expr '(* m (+ (* peri-vxp sin-om) (* peri-vyp cos-om)))}]
      (is (#'core/analytical-genome-valid? ind))
      (is (not (#'core/analytical-genome-valid?
                (assoc ind :qx-expr 'ecc-anom)))))))

(deftest pipeline-only-bootstrap-scores
  (binding [core/*primitive-tier* 5
            core/*primitive-pipeline-only?* true
            core/*both-regimes?* true
            core/*de-driven-search?* true]
    (let [ind (#'core/graded-pipeline-analytical-individual)
          ds (core/scenarios->datasets core/default-scenarios)
          fit (core/calculate-fitness-scenarios ind ds :evaluation :de-driven)]
      (is (#'core/analytical-genome-valid? ind))
      (is (pos? fit)))))

(deftest pipeline-only-bootstrap-uses-pipeline-for-unbound-arm
  (binding [core/*primitive-tier* 5
            core/*primitive-pipeline-only?* true
            core/*both-regimes?* true
            core/*de-driven-search?* true]
    (let [ind (#'core/graded-pipeline-analytical-individual)
          unbound (core/analytical-law-regime-slice ind :unbound)
          hyper (core/scenario-data (last core/default-scenarios))
          metrics (core/evaluate-predictions unbound hyper)]
      (is (= (prims/graded-pipeline-expr :qx-expr) (:qx-expr unbound)))
      (is (= (prims/graded-pipeline-expr :qy-expr) (:qy-expr unbound)))
      (is (< (:mse metrics) 1e-9)))))

(deftest primitive-tier-unlocks-gradually
  (is (= '[e/mean-anomaly] (prims/unlocked-primitive-ops 1)))
  (is (not (some #(= 'e/lab-qx %) (prims/unlocked-primitive-ops 3))))
  (is (some #(= 'e/lab-qx %) (prims/unlocked-primitive-ops 5))))

(deftest compiled-pipeline-expr-matches-oracle
  (binding [core/*primitive-tier* 5
            core/*both-regimes?* false]
    (let [expr (prims/graded-pipeline-expr :qx-expr)
          ind {:strategy :analytical :domain :bound
               :qx-expr expr :qy-expr expr
               :px-expr '(* m peri-vxp) :py-expr '(* m peri-vyp)}
          fns (#'core/compile-analytical-fns ind :cartesian)
          m 1.0 a 1.0 q0x 1.2 q0y 0.0 p0x 0.0 p0y 0.9 t 0.5
          ks (#'core/orbit-state-at-t m a q0x q0y p0x p0y t)]
      (is (prims/expr-primitive-valid? expr 5))
      (is (#'core/analytical-genome-valid? ind))
      (is (some? ks))
      (is (< (Math/abs (- ((:qx fns) t q0x q0y p0x p0y m a) (:qx ks))) 1e-6)))))

(deftest unified-oracle-hyperbola-matches-legacy-slots
  (let [m 1.0 a 1.0 q0x 2.5 q0y 0.3 p0x -0.4 p0y 0.8 t 0.5
        ko (prims/orbit-at-t m a q0x q0y p0x p0y t)
        el (prims/orbital-elements-at-ic m a q0x q0y p0x p0y)
        M  (prims/mean-anomaly (:n-mean el) t (:mean-M0 el))
        F  (prims/anomaly-from-M (:ecc el) (:energy el) M)]
    (is (pos? (:energy el)))
    (is (some? ko))
    (is (< (Math/abs (- F (:hyp-anom ko))) 1e-9))
    (is (< (Math/abs (- (:peri-xp ko)
                        (prims/periapsis-xp (:semi-a el) (:ecc el) (:energy el) F))) 1e-9))))
