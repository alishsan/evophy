(defproject evophy "0.1.0-SNAPSHOT"
  :description "AlphaZero-style learning of physical laws"
  :dependencies [[org.clojure/clojure "1.11.1"]
                 ;; Newer than Emmy’s transitive 0.0.4; avoids abs vs clojure.core warning on 1.11+
                 [org.clojure/math.numeric-tower "0.1.1"]
                 [org.mentat/emmy "0.31.0"]      ;; For the 'Ground Truth' synthetic data
;                 [uncomplicate/neanderthal "0.47.0"] ;; High-performance linear algebra
                 [clj-python/libpython-clj "2.024"]] ;; To bridge to JAX/PyTorch if needed
  :main ^:skip-aot evophy.core
  :target-path "target/%s"
  :profiles {:uberjar {:aot :all}})
