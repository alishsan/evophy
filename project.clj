(defproject evophy "0.1.0-SNAPSHOT"
  :description "Evolve 2D gravitational rate laws and analytical trajectories via GP+MCTS"
  :url "https://github.com/alishsan/evophy"
  :license {:name "Eclipse Public License - v 2.0"
            :url "https://www.eclipse.org/legal/epl-2.0/"}
  :dependencies [[org.clojure/clojure "1.11.1"]
                 ;; Newer than Emmy’s transitive 0.0.4; avoids abs vs clojure.core warning on 1.11+
                 [org.clojure/math.numeric-tower "0.1.1"]
                 [org.mentat/emmy "0.31.0"]      ;; For the 'Ground Truth' synthetic data
;                 [uncomplicate/neanderthal "0.47.0"] ;; High-performance linear algebra
                 [clj-python/libpython-clj "2.024"]] ;; To bridge to JAX/PyTorch if needed
  :main evophy.core
  :target-path "target/%s"
  :aliases {"fitness-progress"  ["run" "-m" "clojure.main" "dev/fitness_progress.clj"]
            "top-expressions"   ["run" "-m" "clojure.main" "dev/top_expressions.clj"]}
  :profiles {:uberjar {:aot [evophy.core]}})
