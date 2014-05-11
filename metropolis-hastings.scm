(define *rejected-samples* 0)
(define *bad-proposals* 0)

(define-structure
  (mh-computation-node
   (constructor make-mh-computation-node
                (operator-instance
                 operator-value
                 running-sampling-logmass
                 sampling-logmass-modifier
                 infeasible?
                 becomes-infeasible?
                 #!optional reused-value)))
  operator-instance
  operator-value
  running-sampling-logmass
  sampling-logmass-modifier
  infeasible?
  becomes-infeasible?
  (reused-value #f))

(define (make-empty-mh-computation-node)
  (make-mh-computation-node
   '() '() '() '() '() '()))

(define-structure
  (mh-state
   (constructor make-mh-state%
                (nsamples burn-in lag output-continuation
                 computation-root computation-path)))
  nsamples
  burn-in
  lag
  output-continuation
  computation-root
  (samples '())
  (samples-drawn 0)
  (samples-recorded 0)
  (common-ancestor)
  ;; mc-value-state and mc-computation-state store the last value
  ;; accepted in the markov chain and its corresponding computation
  ;; path.
  (mc-value-state)
  (mc-computation-state)
  (mc-becomes-infeasible? #f)
  ;; computation-path is a stack of mh-computation-nodes, with the
  ;; special computation-root node at the bottom. Note: this stack
  ;; should /not/ be destructively modified, but rebinding
  ;; computation-path in order to push and pop is fine.
  computation-path
  (becomes-infeasible? #f)
  (named-operator-values (make-weak-eq-hash-table)))

(define (make-mh-state nsamples burn-in lag output-continuation)
  (let ((computation-root
         (make-mh-computation-node
          ;; operator-instance
          #f
          ;; operator-value
          #f
          ;; running-sampling-logmass
          0.0
          ;; sampling-logmass-modifier
          0.0
          ;; infeasible?
          #f
          ;; becomes-infeasible?
          #f)))
    (make-mh-state% nsamples burn-in lag output-continuation
                    computation-root (list computation-root))))

(define (mh-push-computation-node! mh-state node)
  (set-mh-state-computation-path!
   mh-state
   (cons node
         (mh-state-computation-path mh-state))))

(define (instance-resample? instance x)
  (= (instance-mass instance x) 0))

(define (mh-add-to-sampler operator-instance)
  (let* ((name (prob-operator-instance-name operator-instance))
         (default (list 'default))
         (lcrecord (hash-table/get
                    (mh-state-named-operator-values *sampler-state*)
                    name default))
         (resample?
          (cond
           ((null? name) #t) ;; no name means no stored value
           ((eq? lcrecord default) #t) ;; no entry in table
           ;; If the operator instance wasn't part of the previous
           ;; state's computation, then the recorded value is
           ;; stale. Note that mh-resume-with-value-thunk will update
           ;; the record if we don't resample, so we only have to
           ;; check the previous state.
           ((not (mh-in-computation?
                  (local-computation-record-computation-state
                   lcrecord)
                  (mh-state-mc-computation-state *sampler-state*)))
            #t) ;; stored value is from an old branch
           (else (instance-resample?
                  operator-instance
                  (local-computation-record-value lcrecord)))))
         (new-value-thunk
          (if resample?
              (lambda ()
                (instance-sample operator-instance))
              (begin
                ;; (pp "not resampling")
                (let ((new-value
                       (local-computation-record-value lcrecord)))
                  (lambda () new-value))))))
    (mh-resume-with-value-thunk operator-instance new-value-thunk
                                (not resample?))))

(define (mh-resume-with-value-thunk operator-instance value-thunk
                                    #!optional reused?)
  (let ((node-parent
         (car (mh-state-computation-path *sampler-state*)))
        (computation-node (make-empty-mh-computation-node)))
    ;; Push the computation node before evaluating the value thunk,
    ;; because something like mem might want a reference to it
    (mh-push-computation-node! *sampler-state* computation-node)
    (let* ((operator-value (value-thunk))
           (sampling-logmass-modifier
            (instance-logmass operator-instance operator-value))
           (running-sampling-logmass
            (+ (mh-computation-node-running-sampling-logmass node-parent)
               sampling-logmass-modifier))
           (becomes-infeasible?
            (mh-state-becomes-infeasible? *sampler-state*))
           (infeasible?
            (or becomes-infeasible?
                (mh-computation-node-infeasible? node-parent))))
      (set-mh-computation-node-operator-instance!
       computation-node operator-instance)
      (set-mh-computation-node-operator-value!
       computation-node operator-value)
      (set-mh-computation-node-running-sampling-logmass!
       computation-node running-sampling-logmass)
      (set-mh-computation-node-sampling-logmass-modifier!
       computation-node sampling-logmass-modifier)
      (set-mh-computation-node-infeasible?!
       computation-node infeasible?)
      (set-mh-computation-node-becomes-infeasible?!
       computation-node becomes-infeasible?)
      (if (not (default-object? reused?))
          (set-mh-computation-node-reused-value!
           computation-node reused?))
      (set-mh-state-becomes-infeasible?!
       *sampler-state* #f)
      (let ((name (prob-operator-instance-name operator-instance)))
        (if (not (null? name))
            (hash-table/put! (mh-state-named-operator-values *sampler-state*)
                             name
                             (make-local-computation-record operator-value))))
      ((prob-operator-instance-continuation operator-instance)
       operator-value))))

(define (mh-return)
  ((mh-state-output-continuation *sampler-state*)
   ;; samples are collected in reverse order, so put them back before
   ;; returning
   (reverse (mh-state-samples *sampler-state*))))

(define (mh-state-add-sample! state sample)
  (set-mh-state-samples!
   state
   (cons sample
         (mh-state-samples state)))
  (set-mh-state-samples-recorded!
         *sampler-state*
         (+ 1 (mh-state-samples-recorded *sampler-state*))))

(define (mh-maybe-record! sample)
  ;; check burnout, lag, and sample index; add sample if match. Update
  ;; count of samples drawn. Do not resume; handle that elsewhere.
  ;; This also does not update the markov chain's state.
  (if (and (>= (mh-state-samples-drawn *sampler-state*)
               (mh-state-burn-in *sampler-state*))
           (= (modulo (- (mh-state-nsamples *sampler-state*)
                         (mh-state-burn-in *sampler-state*))
                      (mh-state-lag *sampler-state*))))
      (mh-state-add-sample! *sampler-state* sample))
  (set-mh-state-samples-drawn!
   *sampler-state*
   (+ 1 (mh-state-samples-drawn *sampler-state*))))

(define (mh-state-last-sample state)
  (if (null? (mh-state-samples state))
      (error "MH: Tried to get nonexistent last sample")
      (car (mh-state-samples state))))

(define (mh-observe observed)
  (if observed
      'ok
      (set-mh-state-becomes-infeasible?!
       *sampler-state* #t)))

(define (mh-computation-state)
  (mh-state-computation-path *sampler-state*))

(define (mh-in-computation? state computation)
  (let lp ((path computation))
    (cond
     ((eq? path state) #t)
     ((null? path)
      (call-with-current-continuation
       (lambda (k)
         (within-continuation
          (mh-state-output-continuation *sampler-state*)
          (lambda () (k (sampler-in-computation? state)))))))
     (else (lp (cdr path))))))

(define (mh-in-this-computation? state)
  (mh-in-computation?
   state
   (mh-state-computation-path *sampler-state*)))

(define (mh-resample)
  (let ((n-resampling-choices
         (- (length (mh-state-computation-path *sampler-state*)) 1)))
    (cond
     ((= n-resampling-choices 0) ; computation is deterministic
      (error "Can't resample with no probabilistic operators in computation path"))
     (else
      (let* ((resample-index
              (random n-resampling-choices))
             (resample-node
              (list-ref
               (mh-state-computation-path *sampler-state*)
               resample-index))
             (resample-operator
              (mh-computation-node-operator-instance resample-node)))
        (set-mh-state-computation-path!
         *sampler-state*
         (list-tail (mh-state-mc-computation-state *sampler-state*)
                    (+ resample-index 1)))
        (set-mh-state-common-ancestor!
         *sampler-state*
         (car (mh-state-computation-path *sampler-state*)))
        (set-mh-state-becomes-infeasible?!
         *sampler-state* #f)
        (mh-resume-with-value-thunk
         resample-operator
         (lambda ()
           (instance-sample resample-operator))))))))


(define (acceptance-ratio state)
  (cond ((null? (mh-state-samples state))
         ;; This is the initial sample; definitely accept it
         1)
        ((<= (length (mh-state-computation-path state)) 1)
         ;; The computation state doesn't have any probabilistic
         ;; operators in it. It doesn't actually matter whether we
         ;; accept or reject here.
         1)
        ((or (mh-computation-node-infeasible?
              (car (mh-state-mc-computation-state state)))
             (mh-state-mc-becomes-infeasible? state))
         ;; Our previous state was infeasible. Just accept all
         ;; proposals until we find a feasible state. (Note that this
         ;; means that an infeasible sample may be recorded if the
         ;; burn-in period is too low and it is too hard to find an
         ;; initial sample)
         1)
        ((or (mh-computation-node-infeasible?
              (car (mh-state-computation-path state)))
             (mh-state-becomes-infeasible? state))
         ;; Our previous state was feasible, but our new state isn't,
         ;; so the transition probability is 0.
         0)
        (else
         ;; Proceed as normal. The acceptance ratio depends on all the
         ;; nodes for which the proposal distribution differs from the
         ;; distribution to be sampled from, as well as the number of
         ;; nodes in each state's computational history.
         (let ((reused-names (make-eq-hash-table))
               (running-log-ratio 0))
           ;; reused-names will act as a hash set, containing the name
           ;; of every operator for which samples were reused. For the
           ;; current purposes, assume that if a sample is reused
           ;; while proposing x -> x', then that sample would also be
           ;; reused while proposing x' -> x.
           (let lp ((path
                     (mh-state-computation-path state)))
             (if (eq? (car path)
                      (mh-state-common-ancestor state))
                 'done
                 (begin
                   (if (mh-computation-node-reused-value (car path))
                       (let ((name (prob-operator-instance-name
                                    (mh-computation-node-operator-instance
                                     (car path))))
                             (delta-logratio
                              (mh-computation-node-sampling-logmass-modifier
                               (car path))))
                         (hash-table/put! reused-names name #t)
                         (set! running-log-ratio
                               (+ running-log-ratio delta-logratio))))
                   (lp (cdr path)))))
           ;; now that reused-names is populated, find the
           ;; contributions from the original sample.
           (let lp ((path
                     (mh-state-mc-computation-state state)))
             (if (eq? (car path)
                      (mh-state-common-ancestor state))
                 'done
                 (begin
                   (let ((name (prob-operator-instance-name
                                (mh-computation-node-operator-instance
                                 (car path))))
                         (delta-logratio
                          (mh-computation-node-sampling-logmass-modifier
                           (car path))))
                     (if (hash-table/get reused-names name #f)
                         (set! running-log-ratio
                               (- running-log-ratio delta-logratio))))
                   (lp (cdr path)))))
           (let ((length-contribution
                  (/ (- (length (mh-state-mc-computation-state state))
                        1)
                     (- (length (mh-state-computation-path state))
                        1))))
             (* (logmass->mass running-log-ratio)
                length-contribution))))))

(define (mh-query nsamples burn-in lag thunk)
  (if (<= nsamples 0)
      '()
      (call-with-current-continuation
       (lambda (out)
         (fluid-let ((*sampler-state*
                      (make-mh-state nsamples burn-in lag out))
                     (add-to-sampler
                      mh-add-to-sampler)
                     (observe
                      mh-observe)
                     (sampler-computation-state
                      mh-computation-state)
                     (sampler-in-computation?
                      mh-in-this-computation?))
           ;; Draw a sample from the beginning. Then, repeatedly:
           ;; - choose a place to resample from
           ;; - resample from that place
           ;; - compute acceptance ratio
           ;; - accept or reject
           ;; - update state accordingly
           (let lp ((sample (thunk)))
             (let* ((alpha (acceptance-ratio *sampler-state*))
                    ;; note: (random 1.0) returns a number between
                    ;; zero inclusive and one exclusive. if alpha is 0
                    ;; and we draw a 0, we need to make sure to
                    ;; reject. So we accept if the acceptance ratio is
                    ;; strictly greater than the number drawn.
                    (accepted (< (random 1.0) alpha)))
               (if (= alpha 0)
                   (set! *bad-proposals* (+ *bad-proposals* 1)))
               (if accepted
                   (begin
                     (mh-maybe-record! sample)
                     (set-mh-state-mc-value-state! *sampler-state* sample)
                     (set-mh-state-mc-computation-state!
                      *sampler-state*
                      (mh-state-computation-path *sampler-state*))
                     (set-mh-state-mc-becomes-infeasible?!
                      *sampler-state*
                      (mh-state-becomes-infeasible? *sampler-state*)))
                   (begin (mh-maybe-record! (mh-state-mc-value-state
                                             *sampler-state*))
                          (set! *rejected-samples*
                                (+ *rejected-samples* 1))))
               (cond
                ((>= (mh-state-samples-recorded *sampler-state*)
                     (mh-state-nsamples *sampler-state*))
                 (mh-return)) ; we're finished
                ((<=
                  (length (mh-state-computation-path *sampler-state*))
                  1)
                 ;; no probabilistic operators, so no need to resample
                 (lp sample)) 
                (else (mh-resample))))))))))