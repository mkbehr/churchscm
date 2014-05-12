(define-structure
  (rejection-state
   (constructor make-rejection-state
                (nsamples cutoff output-continuation)))
  (top-level '())
  (samples '())
  nsamples
  cutoff
  (computation-state (generate-uninterned-symbol))
  output-continuation)

(define (rejection-add-to-sampler operator-instance)
  (if (null? (rejection-state-top-level *sampler-state*))
      (set-rejection-state-top-level!
       *sampler-state*
       operator-instance))
  (let ((value (instance-sample operator-instance)))
    (if (and (prob-operator-instance-constrained? operator-instance)
             (< (instance-mass instance value)
               (random 1.0)))
        (rejection-reject)
        ((prob-operator-instance-continuation operator-instance)
         value))))

(define (rejection-return)
  ((rejection-state-output-continuation *sampler-state*)
   ;; samples are collected in reverse order, so put them back before
   ;; returning
   (reverse (rejection-state-samples *sampler-state*))))

(define (rejection-state-add-sample! state sample)
  (set-rejection-state-samples!
   state
   (cons sample
         (rejection-state-samples state)))
  (set-rejection-state-nsamples!
   state
   (- (rejection-state-nsamples state) 1))
  (set-rejection-state-cutoff!
   state
   (- (rejection-state-cutoff state) 1)))

(define (rejection-reject)
  (set-rejection-state-cutoff!
   *sampler-state*
   (- (rejection-state-cutoff *sampler-state*) 1))
  (if (<= ;; have we tried too many times?
       (rejection-state-cutoff *sampler-state*)
       0)
      (rejection-return)
      (rejection-restart)))

(define (rejection-observe observed)
  (if observed
      'ok
      (rejection-reject)))

(define (rejection-computation-state)
  (rejection-state-computation-state *sampler-state*))

(define (rejection-in-computation? state)
  (if (eq? state (sampler-computation-state))
      #t
      ;; If e.g. a mem table entry wasn't defined in our current
      ;; branch, we may still be doing nested queries, so check to see
      ;; if it would have been in our current branch on the next
      ;; sampling level up. The output continuation should contain all
      ;; the information we need in order to do that.
      (call-with-current-continuation
       (lambda (k)
         (within-continuation
          (rejection-state-output-continuation *sampler-state*)
          (lambda () (k (sampler-in-computation? state))))))))

(define (rejection-restart)
  (let ((top-level (rejection-state-top-level *sampler-state*)))
    (set-rejection-state-computation-state!
     *sampler-state* (generate-uninterned-symbol))
    ((prob-operator-instance-continuation top-level)
     (instance-sample top-level))))

(define (rejection-query samples cutoff thunk)
  (if (or (<= samples 0) (<= cutoff 0))
      '()
      (call-with-current-continuation
         (lambda (out)
           (fluid-let ((*sampler-state*
                        (make-rejection-state samples cutoff out))
                       (add-to-sampler
                        rejection-add-to-sampler)
                       (observe
                        rejection-observe)
                       (sampler-computation-state
                        rejection-computation-state)
                       (sampler-in-computation?
                        rejection-in-computation?))
             (let lp ((sample (thunk)))
               (rejection-state-add-sample! *sampler-state* sample)
               (if (or
                    (<= (rejection-state-nsamples *sampler-state*) 0)
                    (<= (rejection-state-cutoff *sampler-state*) 0))
                   (rejection-return)
                   (if (null?
                        (rejection-state-top-level *sampler-state*))
                       (lp (thunk))
                       (rejection-restart)))))))))