;;; Tests for system without any inference running

(flip)
;; expected output: 1 or 0
(observe #t)
;; expected output: not an error
(observe #f)
;; expected output: error

;;; Continuous distributions

(cont-uniform 0 10)
;; Expected value: inexact number between 0 and 10
(pspec-logmass cont-uniform-spec (list 0 10) 8)
(pspec-logmass cont-uniform-spec (list 0 10) 9)
;; Expected value: the same small value
(pspec-logmass cont-uniform-spec (list 0 10) 20)
;; Expected value: -inf
(normal 0 1)
;; Expected value: sample from the standard normal distribution
(pspec-logmass normal-spec (list 0 1) 0)
;; Expected value: something small
(pspec-logmass normal-spec (list 0 1) 100)
;; Expected value: something extremely small (this will be -inf
;; barring a very precise computer)

;;; Tests for rejection sampling

(rejection-query
 400 4000
 (lambda ()
   (pp "This should only print once")
   (flip)))
;; expected output: print "This should only print once", return list
;; of 400 1s and 0s, evenly distributed

(equal-histogram
 (rejection-query
  400 4000
  (lambda ()
    (flip))))

(rejection-query
 10 100
 (lambda () 1))
;; expected output: list of 10 1

;;; Tests for rejection sampling with observation

(rejection-query
 20 200
 (lambda ()
   (pp "This should only print once")
   (let ((result1 (flip))
         (result2 (flip)))
     (observe (= (+ result1 result2) 1))
     (list result1 result2))))
;; expected output: print "This should only print once", return list
;; of 20 (1 0) and (0 1), evenly distributed

(rejection-query
 10 100
 (lambda ()
   (let ((result (flip)))
     (observe #f)
     result)))
;; expected output: empty list

(rejection-query
 20 200
 (lambda ()
   (let ((result1 (bernoulli 2/3))
         (result2 (bernoulli 1/3)))
     (observe (= (+ result1 result2) 1))
     (list result1 result2))))
;; expected output: list of 20 (1 0) and (0 1), with a ~4:1 ratio

(let ((n 1000))
  (/ (fold-left
      + 0
      (rejection-query
       n n
       (lambda () (normal 0 1))))
     n))
;; Mean of n standard normal variables. Expected output: a small
;; number.

;;; Tests for rejection sampling with memoization

(rejection-query
 20 200
 (lambda ()
   (let ((f (mem flip)))
     (list (f) (f)))))
;; expected output: list of 20 (1 1) and (0 0), ~evenly distributed

(let ((f (mem flip)))
  (rejection-query
   20 200
   (lambda ()
     (list (f) (f)))))
;; expected output: list of 20 (1 1) and (0 0), ~evenly distributed

(let ((f (mem flip)))
  (f)
  (rejection-query
   20 200
   (lambda ()
     (list (f) (f)))))
;; expected output: list of 20 (1 1) or 20 (0 0)

;;; Test for nested rejection sampling and mem

(pp
 (rejection-query
  5 200
  (lambda ()
    (let* ((f (mem bernoulli))
           (f-result (f 49/100)))
      (rejection-query
       5 200
       (lambda ()
         (list
          f-result (f 49/100)
          (f 51/100) (f 51/100))))))))
;; expected output: list of five lists of five lists of four numbers:
;; (((w x y z))), each 1 or 0
;;
;; in each bottom-level list, w = x and y = z
;;
;; in each middle level list, every w should be the same as every
;; other w, but y should be approximately evenly distributed
;;
;; in the top-level list, w should be approximately evenly distributed
;; (actually 51% 0 for w and 51% 1 for y, but that won't be visible on
;; this level)


;;; Tests for metropolis-hastings

(let ((f (mem flip)))
  (mh-query
   20 10 10
   (lambda ()
     (list (f) (f)))))

(define (simple-bayes-net observation)
  (let* ((cloudy (bernoulli 0.5))
         (sprinkler-prob
          (if (= cloudy 1)
              0.1
              0.5))
         (sprinkler (bernoulli sprinkler-prob))
         (rain-prob
          (if (= cloudy 1)
              0.8
              0.2))
         (rain (bernoulli rain-prob))
         (wet-grass-prob
          (cond
           ((and (= sprinkler 1) (= rain 1))
            0.99)
           ((and (= sprinkler 1) (= rain 0))
            0.9)
           ((and (= sprinkler 0) (= rain 1))
            0.9)
           (else
            0.01)))
         (wet-grass (bernoulli wet-grass-prob)))
    (if observation
        (for-each
         (lambda (val ob)
           (if ob
               (observe (= val ob))))
         (list cloudy sprinkler rain wet-grass)
         observation))
    (list cloudy sprinkler rain wet-grass)))

(equal-histogram
 (rejection-query
  10000 100000
  (lambda ()
    (simple-bayes-net '(#f #f #f 1)))))

(equal-histogram
 (mh-query
  10000 1000 10
  (lambda ()
    (simple-bayes-net '(#f #f #f 1)))))
;; these two should return similar values



(with-timings
 (lambda ()
   (equal-histogram
    (rejection-query
     10000 100000
     (lambda ()
       (simple-bayes-net '(#f #f #f 1))))))
 (lambda (run-time gc-time real-time)
   (write (internal-time/ticks->seconds run-time))
   (write-char #\space)
   (write (internal-time/ticks->seconds gc-time))
   (write-char #\space)
   (write (internal-time/ticks->seconds real-time))
   (newline)))

(with-timings
 (lambda ()
   (equal-histogram
    (mh-query
     10000 100 10
     (lambda ()
       (simple-bayes-net '(#f #f #f 1))))))
 (lambda (run-time gc-time real-time)
   (write (internal-time/ticks->seconds run-time))
   (write-char #\space)
   (write (internal-time/ticks->seconds gc-time))
   (write-char #\space)
   (write (internal-time/ticks->seconds real-time))
   (newline)))
;; MCMC is hopefully faster! Though note that the speed of MCMC is
;; determined entirely by the sampler's parameters and the time it
;; takes to draw a sample, while the speed of rejection sampling also
;; depends on the probability of rejecting a sample.

(let ((n 1000))
  (/ (fold-left
      + 0
      (mh-query
       n 100 10
       (lambda () (normal 0 1))))
     n))
;; Mean of standard normal variables again. Expected output: a small
;; number.

(define (named-bayes-net observation)
  (let* ((cloudy ((named-operator bernoulli 'cloudy) 0.5))
         (sprinkler-prob
          (if (= cloudy 1)
              0.1
              0.5))
         (sprinkler ((named-operator bernoulli 'sprinkler)
                     sprinkler-prob))
         (rain-prob
          (if (= cloudy 1)
              0.8
              0.2))
         (rain ((named-operator bernoulli 'rain)
                rain-prob))
         (wet-grass-prob
          (cond
           ((and (= sprinkler 1) (= rain 1))
            0.99)
           ((and (= sprinkler 1) (= rain 0))
            0.9)
           ((and (= sprinkler 0) (= rain 1))
            0.9)
           (else
            0.01)))
         (wet-grass ((named-operator bernoulli 'wet-grass)
                     wet-grass-prob)))
    (if observation
        (for-each
         (lambda (val ob)
           (if ob
               (observe (= val ob))))
         (list cloudy sprinkler rain wet-grass)
         observation))
    (list cloudy sprinkler rain wet-grass)))

(equal-histogram
 (rejection-query
  10000 100000
  (lambda ()
    (simple-bayes-net '(#f #f 1 1)))))

(equal-histogram
 (rejection-query
  10000 100000
  (lambda ()
    (named-bayes-net '(#f #f 1 1)))))

(define *rejected-samples* 0)
(define *bad-proposals* 0)
(equal-histogram
 (mh-query
  10000 10000 100
  (lambda ()
    (simple-bayes-net '(#f #f 1 1)))))
*rejected-samples*
*bad-proposals*

(define *rejected-samples* 0)
(define *bad-proposals* 0)
(equal-histogram
 (mh-query
  10000 10000 100
  (lambda ()
    (named-bayes-net '(#f #f 1 1)))))
*rejected-samples*
*bad-proposals*
;; Sampling from named-bayes-net should give fewer rejected proposals
;; and bad proposals than sampling from simple-bayes-net. Also,
;; simple-bayes-net should give equal numbers of rejected proposals
;; and bad proposals, but named-bayes-net should reject some feasible
;; proposals.

(define (constrained-bayes-net observation)
  (let ((cloudy-constraint (list-ref observation 0))
        (sprinkler-constraint (list-ref observation 1))
        (rain-constraint (list-ref observation 2))
        (wet-grass-constraint (list-ref observation 3)))
    (let* ((cloudy-operator
            (if cloudy-constraint
                (constrain bernoulli cloudy-constraint)
                (named-operator bernoulli 'cloudy)))
           (cloudy (cloudy-operator 0.5))
           (sprinkler-operator
            (if sprinkler-constraint
                (constrain bernoulli sprinkler-constraint)
                (named-operator bernoulli 'sprinkler)))
           (sprinkler-prob
            (if (= cloudy 1)
                0.1
                0.5))
           (sprinkler (sprinkler-operator
                       sprinkler-prob))
           (rain-operator
            (if rain-constraint
                (constrain bernoulli rain-constraint)
                (named-operator bernoulli 'rain)))
           (rain-prob
            (if (= cloudy 1)
                0.8
                0.2))
           (rain (rain-operator
                  rain-prob))
           (wet-grass-operator
            (if wet-grass-constraint
                (constrain bernoulli wet-grass-constraint)
                (named-operator bernoulli 'wet-grass)))
           (wet-grass-prob
            (cond
             ((and (= sprinkler 1) (= rain 1))
              0.99)
             ((and (= sprinkler 1) (= rain 0))
              0.9)
             ((and (= sprinkler 0) (= rain 1))
              0.9)
             (else
              0.01)))
           (wet-grass (wet-grass-operator
                       wet-grass-prob)))
      (list cloudy sprinkler rain wet-grass))))

(equal-histogram
 (rejection-query
  10000 100000
  (lambda ()
    (constrained-bayes-net '(#f #f 1 1)))))

(define *rejected-samples* 0)
(define *bad-proposals* 0)
(equal-histogram
 (mh-query
  10000 10000 10
  (lambda ()
    (constrained-bayes-net '(#f #f 1 1)))))
*rejected-samples*
*bad-proposals*
;; This should tend to have fewer rejected samples than the named
;; bayes net, and it should have no bad proposals at all.