(ns evophy.primitives-test
  (:require [clojure.test :refer :all]
            [evophy.core :as core]
            [evophy.primitives :as prims]))

(deftest mean-anomaly-and-anomaly-match-kepler-u
  (let [m 1.0 a 1.0 q0x 1.2 q0y 0.0 p0x 0.0 p0y 0.9 t 0.5
        el (prims/orbital-elements-at-ic m a q0x q0y p0x p0y)
        M  (prims/mean-anomaly (:n-mean el) t (:kepler-M0 el))
        u  (prims/anomaly-from-M (:ecc el) (:energy el) M)
        ko (#'core/kepler-orbit-at-t m a q0x q0y p0x p0y t)]
    (is (some? el))
    (is (some? ko))
    (is (< (Math/abs (- u (:kepler-u ko))) 1e-9))))

(deftest periapsis-pipeline-matches-kepler-position
  (let [m 1.0 a 1.0 q0x 1.2 q0y 0.0 p0x 0.0 p0y 0.9 t 0.5
        el (prims/orbital-elements-at-ic m a q0x q0y p0x p0y)
        M  (prims/mean-anomaly (:n-mean el) t (:kepler-M0 el))
        u  (prims/anomaly-from-M (:ecc el) (:energy el) M)
        xp (prims/periapsis-xp (:semi-a el) (:ecc el) (:energy el) u)
        yp (prims/periapsis-yp (:semi-b el) (:bh el) (:energy el) u)
        ko (#'core/kepler-orbit-at-t m a q0x q0y p0x p0y t)]
    (is (< (Math/abs (- xp (:kepler-xp ko))) 1e-9))
    (is (< (Math/abs (- yp (:kepler-yp ko))) 1e-9))))

(deftest graded-pipeline-depth-capped
  (let [deep '(e/lab-qx cos-om sin-om
              (e/periapsis-xp semi-a ecc energy
                (e/anomaly ecc energy
                  (e/mean-anomaly n-mean t kepler-M0)))
              (e/periapsis-yp semi-b bh energy kepler-u))]
    (is (prims/expr-primitive-valid? deep 5))
    (is (not (prims/expr-primitive-valid? deep 2)))))

(deftest primitive-tier-unlocks-gradually
  (is (= '[e/mean-anomaly] (prims/unlocked-primitive-ops 1)))
  (is (not (some #(= 'e/lab-qx %) (prims/unlocked-primitive-ops 3))))
  (is (some #(= 'e/lab-qx %) (prims/unlocked-primitive-ops 5))))

(deftest compiled-pipeline-expr-matches-oracle
  (binding [core/*primitive-tier* 5
            core/*both-regimes?* false]
    (let [expr (prims/kepler-pipeline-expr :qx-expr)
          ind {:strategy :analytical :domain :bound
               :qx-expr expr :qy-expr expr
               :px-expr '(* m kepler-vxp) :py-expr '(* m kepler-vyp)}
          fns (#'core/compile-analytical-fns ind :cartesian)
          m 1.0 a 1.0 q0x 1.2 q0y 0.0 p0x 0.0 p0y 0.9 t 0.5
          ks (#'core/kepler-state-at-t m a q0x q0y p0x p0y t)]
      (is (prims/expr-primitive-valid? expr 5))
      (is (#'core/analytical-genome-valid? ind))
      (is (some? ks))
      (is (< (Math/abs (- ((:qx fns) t q0x q0y p0x p0y m a) (:qx ks))) 1e-6)))))

(deftest unified-oracle-hyperbola-matches-legacy-slots
  (let [m 1.0 a 1.0 q0x 2.5 q0y 0.3 p0x -0.4 p0y 0.8 t 0.5
        ko (prims/kepler-orbit-at-t m a q0x q0y p0x p0y t)
        el (prims/orbital-elements-at-ic m a q0x q0y p0x p0y)
        M  (prims/mean-anomaly (:n-mean el) t (:kepler-M0 el))
        F  (prims/anomaly-from-M (:ecc el) (:energy el) M)]
    (is (pos? (:energy el)))
    (is (some? ko))
    (is (< (Math/abs (- F (:kepler-F ko))) 1e-9))
    (is (< (Math/abs (- (:kepler-xp ko)
                        (prims/periapsis-xp (:semi-a el) (:ecc el) (:energy el) F))) 1e-9))))
