(ns evophy.symbolic
  "Bridge between evophy GP S-expressions and Emmy symbolic expressions.
   Used to validate and simplify catalog blocks at load time (not per fitness eval)."
  (:require [clojure.walk :as walk]
            [emmy.abstract.number :as an]
            [emmy.expression :as ex]
            [emmy.generic :as g]))

(def ^:private emmy-sym (comp an/literal-number symbol))

(defn- evophy-op->emmy [op]
  (cond
    (= op '+) g/+
    (= op '-) g/-
    (= op '*) g/*
    (= op 'e/square) g/square
    (= op 'e/sin) g/sin
    (= op 'e/cos) g/cos
    (= op 'e/sqrt) g/sqrt
    (or (= op 'e/div) (= op 'e//)) g//
    :else nil))

(defn evophy-expr->emmy
  "Convert an evophy expression tree to an Emmy symbolic form."
  [expr]
  (cond
    (number? expr) (double expr)
    (symbol? expr) (emmy-sym expr)
    (and (sequential? expr) (seq expr))
    (let [[op & args] expr
          emmy-op (evophy-op->emmy op)]
      (if emmy-op
        (apply emmy-op (map evophy-expr->emmy args))
        (throw (ex-info "Unsupported evophy op for Emmy" {:op op}))))
    :else expr))

(defn- emmy-op->evophy-sym [x]
  (cond
    (instance? clojure.lang.Named x) (symbol x)
    :else x))

(defn- emmy-expr->evophy*
  [x]
  (cond
    (number? x) x
    (ex/literal? x) (emmy-op->evophy-sym (ex/expression-of x))
    (and (sequential? x) (seq x))
    (let [[op & args] x]
      (case op
        g/+ (cons '+ (map emmy-expr->evophy* args))
        g/- (if (= 1 (count args))
              (list '- (emmy-expr->evophy* (first args)))
              (cons '- (map emmy-expr->evophy* args)))
        g/* (cons '* (map emmy-expr->evophy* args))
        g/square (list 'e/square (emmy-expr->evophy* (first args)))
        g/sin (list 'e/sin (emmy-expr->evophy* (first args)))
        g/cos (list 'e/cos (emmy-expr->evophy* (first args)))
        g/sqrt (list 'e/sqrt (emmy-expr->evophy* (first args)))
        g// (list 'e/div (emmy-expr->evophy* (first args)) (emmy-expr->evophy* (second args)))
        (throw (ex-info "Unsupported Emmy op in round-trip" {:op op}))))
    :else x))

(defn emmy-expr->evophy [x]
  (walk/postwalk
   (fn [node]
     (if (ex/literal? node)
       (emmy-op->evophy-sym (ex/expression-of node))
       node))
   (emmy-expr->evophy* x)))

(defn simplify-evophy-expr
  "Best-effort algebraic simplification via Emmy; returns original on failure."
  [expr]
  (try
    (-> expr evophy-expr->emmy g/simplify emmy-expr->evophy)
    (catch Throwable _ expr)))

(defn symbolic-validate-expr?
  "True when expr converts to Emmy (e/if arms validated recursively)."
  [expr]
  (try
    (cond
      (or (number? expr) (symbol? expr)) true
      (and (sequential? expr) (= 'e/if (first expr)))
      (and (= 4 (count expr))
           (= '(neg? energy) (second expr))
           (symbolic-validate-expr? (nth expr 2))
           (symbolic-validate-expr? (nth expr 3)))
      (and (sequential? expr) (seq expr))
      (boolean (evophy-expr->emmy expr))
      :else false)
    (catch Throwable _ false)))
