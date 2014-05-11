(define *pi* (acos -1))
(define +inf (/ 1.0 0.0))
(define -inf (/ -1.0 0.0))
(define NaN (+ +inf -inf))

(define (equal-histogram samples)
  (define hist-table (make-equal-hash-table))
  (define not-in-table (list 'not-in-table))
  (for-each
   (lambda (sample)
     (if (eq? (hash-table/get hist-table sample not-in-table)
              not-in-table)
         (hash-table/put! hist-table sample 1)
         (hash-table/put!
          hist-table sample
          (+ 1 (hash-table/get hist-table sample not-in-table)))))
   samples)
  (hash-table->alist hist-table))