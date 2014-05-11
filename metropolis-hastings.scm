(define *rejected-samples* 0)
(define *bad-proposals* 0)

(define-structure
  (mh-computation-node
   (constructor make-mh-computation-node
                (operator-instance
                 operator-value
                 running-sampling-logmass
                 running-observed-logmass
                 sampling-logmass-modifier
                 observed-logmass-modifier
                 #!optional reused-value)))
  operator-instance
  operator-value
  running-sampling-logmass
  ;; note: observed mass is not actually a proper probability mass,
  ;; because it is less than or equal to sampling mass.
  running-observed-logmass
  sampling-logmass-modifier
  observed-logmass-modifier
  (reused-value #f))

(define (make-empty-mh-computation-node)
  (make-mh-computation-node
   '() '() '() '() '() '()))

(define-structure
  (mh-state
   (constructor make-mh-state%
                ;; TODO figure out what lag is canonically called -
                ;; the number of samples to throw out between every
                ;; sample we draw
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
  ;; mc-value-state and mc-computation-state store the last value
  ;; accepted in the markov chain and its corresponding computation
  ;; path.
  (mc-value-state)
  (mc-computation-state)
  ;; There's currently some gymnastics required to keep track of this
  ;; and reset running-observed-mass-modifier at appropriate times,
  ;; because any observations made after the last probabilistic
  ;; operator are not reflected in the computation state structure.
  (mc-observed-logmass-modifier 0)
  ;; computation-path is a stack of mh-computation-nodes, with the
  ;; special computation-root node at the bottom. Note: this stack
  ;; should /not/ be destructively modified, but rebinding
  ;; computation-path in order to push and pop is fine.
  computation-path
  (running-observed-logmass-modifier 0)
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
          ;; running-observed-logmass
          0.0
          ;; sampling-logmass-modifier
          0.0
          ;; observed-logmass-modifier
          0.0)))
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
           ;; TODO explain this if it works
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
           (observed-logmass-modifier
            (+ sampling-logmass-modifier
               (mh-state-running-observed-logmass-modifier
                *sampler-state*)))
           (running-sampling-logmass
            (+ (mh-computation-node-running-sampling-logmass node-parent)
               sampling-logmass-modifier))
           (running-observed-logmass
            (+ (mh-computation-node-running-observed-logmass node-parent)
               observed-logmass-modifier)))
      (set-mh-computation-node-operator-instance!
       computation-node operator-instance)
      (set-mh-computation-node-operator-value!
       computation-node operator-value)
      (set-mh-computation-node-running-sampling-logmass!
       computation-node running-sampling-logmass)
      (set-mh-computation-node-running-observed-logmass!
       computation-node running-observed-logmass)
      (set-mh-computation-node-sampling-logmass-modifier!
       computation-node sampling-logmass-modifier)
      (set-mh-computation-node-observed-logmass-modifier!
       computation-node observed-logmass-modifier)
      (if (not (default-object? reused?))
          (set-mh-computation-node-reused-value!
           computation-node reused?))
      (set-mh-state-running-observed-logmass-modifier!
       *sampler-state* 0)
      (let ((name (prob-operator-instance-name operator-instance)))
        (if (not (null? name))
            (hash-table/put! (mh-state-named-operator-values *sampler-state*)
                             name
                             (make-local-computation-record operator-value))))
      ((prob-operator-instance-continuation operator-instance)
       operator-value))))

;; TODO many functions here look the same as in rejection
;; sampling. Consider abstracting. Will mark with "TODO duplicated
;; code; see above".
(define (mh-return)
  ((mh-state-output-continuation *sampler-state*)
   ;; samples are collected in reverse order, so put them back before
   ;; returning
   (reverse (mh-state-samples *sampler-state*))))

;; TODO duplicated code; see above
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
  (set-mh-state-running-observed-logmass-modifier!
   *sampler-state*
   (+ (mh-state-running-observed-logmass-modifier *sampler-state*)
      (mass->logmass observed))))

(define (mh-computation-state)
  (mh-state-computation-path *sampler-state*))

(define (mh-in-computation? state computation)
  (let lp ((path computation))
    (cond
     ((eq? path state) #t)
     ((null? path)
      ;; TODO duplicated code; see above
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
         (list-tail (mh-state-computation-path *sampler-state*)
                    (+ resample-index 1)))
        (set-mh-state-running-observed-logmass-modifier!
         *sampler-state* 0)
        (mh-resume-with-value-thunk
         resample-operator
         (lambda ()
           (instance-sample resample-operator))))))))

;; TODO this could use a better name
(define (log-sampling-mass-difference old-path new-path)
  (define out 0)
  ;; Note: for this purpose, we assume that whether or not an operator
  ;; is resampled is deterministic and symmetric.
  (define old-path-entries (make-eq-hash-table))
  (for-each
   (lambda (node)
     (hash-table/put! old-path-entries node #t))
   old-path)
  ;; Assume the common ancestor exists. It will unless something has
  ;; gone wrong.
  (let ((common-ancestor
         (let lp ((path new-path))
           (if (hash-table/get old-path-entries (car path) #f)
               (car path)
               (lp (cdr path))))))
    (define name-table (make-eq-hash-table))
    ;; Go through all the nodes in the old computation up to the
    ;; common ancestor: if they don't have a name, add them to the
    ;; bias right away. Otherwise, store them in the name table. As
    ;; always, assume that no two operators have the same name.
    (let lp ((path old-path))
      (if (eq? (car path) common-ancestor)
          'done
          (let ((name (prob-operator-instance-name
                       (mh-computation-node-operator-instance
                        (car path)))))
            (if (null? name)
                (set! out
                      (+ out
                         (mh-computation-node-sampling-logmass-modifier
                          (car path))))
                (hash-table/put! name-table name (car path)))
            (lp (cdr path)))))
    ;; Now go through the nodes in the new computation up to the
    ;; common ancestor. Check to see if they have a counterpart in the
    ;; name table. If the counterpart exists and we didn't resample,
    ;; then by our earlier assumption that resampling is deterministic
    ;; and symmetric, we can ignore them both. Otherwise, proceed as
    ;; normal: subtract them from the bias and leave the entry in the
    ;; name table to be taken into account at the next step.
    (let lp ((path new-path))
      (if (eq? (car path) common-ancestor)
          'done
          (let* ((name (prob-operator-instance-name
                        (mh-computation-node-operator-instance
                         (car path))))
                 (lookup-object
                  (hash-table/get name-table name #f)))
            (if (and lookup-object
                     (mh-computation-node-reused-value
                      (car path)))
                (hash-table/remove! name-table name)
                (set! out
                      (- out
                         (mh-computation-node-sampling-logmass-modifier
                          (car path))))))))
    ;; Finally, take into account any remaining nodes in the name
    ;; table.
    (for-each
     (lambda (node)
       (set! out
             (+ out
                (mh-computation-node-sampling-logmass-modifier
                 node))))
     (hash-table/datum-list name-table)))
  out)


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
             (let* ((acceptance-ratio
                     (cond
                      ((null? (mh-state-samples *sampler-state*))
                       ;; This is the initial sample; definitely accept it
                       1)
                      ((<= (length
                            (mh-state-computation-path
                             *sampler-state*))
                           1)
                       ;; The computation state doesn't have any
                       ;; probabilistic operators in it. Calculating
                       ;; the acceptance ratio normally would divide
                       ;; by zero, but all samples are the same, so we
                       ;; can just accept.
                       1)
                      ((= (+
                           (mh-computation-node-running-observed-logmass
                            (car (mh-state-mc-computation-state
                                  *sampler-state*)))
                           (mh-state-mc-observed-logmass-modifier
                            *sampler-state*))
                          -inf)
                       ;; Our previous state was infeasible. This
                       ;; would again result in division by zero, so
                       ;; just accept all proposals until we find a
                       ;; feasible state. (Note that this means that
                       ;; an infeasible sample may be recorded if the
                       ;; burn-in period is too low and it is too hard
                       ;; to find an initial sample)
                       1)
                      (else
                       ;; Proceed as normal. The acceptance ratio is
                       ;; (P(x')/P(x)) * (g(x' -> x) / g(x -> x'))
                       ;; (capped at 1, but this happens implicitly).
                       (let ((logmass-difference
                              ;; Note: this will divide by zero if a
                              ;; sample with sampling mass of 0 is
                              ;; ever accepted. This should never
                              ;; happen, but it is important to take
                              ;; care that the initial sample is
                              ;; feasible.
                              (-
                               ;; get the running mass we just
                               ;; computed, making sure to factor in
                               ;; any last-minute observations...
                               (+
                                (mh-computation-node-running-observed-logmass
                                 (car (mh-state-computation-path
                                       *sampler-state*)))
                                (mh-state-running-observed-logmass-modifier
                                 *sampler-state*))
                               ;; and do the same for our previous
                               ;; state.
                               (+
                                (mh-computation-node-running-observed-logmass
                                 (car (mh-state-mc-computation-state
                                       *sampler-state*)))
                                (mh-state-mc-observed-logmass-modifier
                                 *sampler-state*))))
                             ;; Assuming that we make a uniform
                             ;; choice of which operator to resample
                             ;; from, the bias ratio depends only on
                             ;; the sampling mass and the length of
                             ;; each state's computation path. See my
                             ;; notes from may 3 (TODO: see the
                             ;; writeup) for why.
                             (logbias-difference
                              (+
                               (mass->logmass
                                (/ (- (length
                                       (mh-state-mc-computation-state
                                        *sampler-state*))
                                      1)
                                   (- (length
                                       (mh-state-computation-path
                                        *sampler-state*))
                                      1)))
                               (log-sampling-mass-difference
                                (mh-state-mc-computation-state
                                 *sampler-state*)
                                (mh-state-computation-path
                                 *sampler-state*)))))
                         ;; DEBUG
                         ;; (pp (list logmass-difference logbias-difference))
                         (logmass->mass
                          (+ logmass-difference logbias-difference))))))
                    ;; note: (random 1.0) returns a number between
                    ;; zero inclusive and one exclusive. if alpha is 0
                    ;; and we draw a 0, we need to make sure to
                    ;; reject. So we accept if the acceptance ratio is
                    ;; strictly greater than the number drawn.
                    (accepted (< (random 1.0) acceptance-ratio)))
               (if (= acceptance-ratio 0)
                   (set! *bad-proposals* (+ *bad-proposals* 1)))
               ;; DEBUG
               (if #f
                   (pp (list (mh-state-mc-value-state *sampler-state*)
                             sample
                             acceptance-ratio
                             (if (and
                                  (mh-state-mc-computation-state
                                   *sampler-state*)
                                  (mh-state-computation-path
                                   *sampler-state*))
                                 (-
                                  (+
                                   (mh-computation-node-running-observed-logmass
                                    (car (mh-state-computation-path
                                          *sampler-state*)))
                                   (mh-state-running-observed-logmass-modifier
                                    *sampler-state*))
                                  (mh-computation-node-running-sampling-logmass
                                   (car (mh-state-mc-computation-state
                                         *sampler-state*)))))
                             (if (mh-state-mc-computation-state
                                  *sampler-state*)
                                 (exp (mh-computation-node-running-sampling-logmass
                                       (car (mh-state-mc-computation-state
                                             *sampler-state*)))))
                             (if (mh-state-computation-path
                                  *sampler-state*)
                                 (exp (mh-computation-node-running-sampling-logmass
                                       (car (mh-state-computation-path
                                             *sampler-state*))))))))
               (if accepted
                   (begin
                     (mh-maybe-record! sample)
                     (set-mh-state-mc-value-state! *sampler-state* sample)
                     (set-mh-state-mc-computation-state! *sampler-state*
                                                         (mh-state-computation-path *sampler-state*))
                     (set-mh-state-mc-observed-logmass-modifier!
                      *sampler-state*
                      (mh-state-running-observed-logmass-modifier *sampler-state*)))
                   (begin (mh-maybe-record! (mh-state-mc-value-state
                                             *sampler-state*))
                          ;; DEBUG?
                          (set! *rejected-samples* (+
                                                    *rejected-samples*
                                                    1))))
               (cond
                ((>= (mh-state-samples-recorded *sampler-state*)
                     (mh-state-nsamples *sampler-state*))
                 (mh-return)) ; we're finished
                ((<= (length (mh-state-computation-path *sampler-state*)) 1)
                 (lp sample)) ; no probabilistic operators, so no
                                        ; need to resample
                (else (mh-resample))))))))))