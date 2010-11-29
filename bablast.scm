(define BLOSUM62 (list->vector '(
  4  -1  -2  -2   0  -1  -1   0  -2  -1  -1  -1  -1  -2  -1   1   0  -3  -2   0 
 -1   5   0  -2  -3   1   0  -2   0  -3  -2   2  -1  -3  -2  -1  -1  -3  -2  -3 
 -2   0   6   1  -3   0   0   0   1  -3  -3   0  -2  -3  -2   1   0  -4  -2  -3 
 -2  -2   1   6  -3   0   2  -1  -1  -3  -4  -1  -3  -3  -1   0  -1  -4  -3  -3 
  0  -3  -3  -3   9  -3  -4  -3  -3  -1  -1  -3  -1  -2  -3  -1  -1  -2  -2  -1 
 -1   1   0   0  -3   5   2  -2   0  -3  -2   1   0  -3  -1   0  -1  -2  -1  -2 
 -1   0   0   2  -4   2   5  -2   0  -3  -3   1  -2  -3  -1   0  -1  -3  -2  -2 
  0  -2   0  -1  -3  -2  -2   6  -2  -4  -4  -2  -3  -3  -2   0  -2  -2  -3  -3 
 -2   0   1  -1  -3   0   0  -2   8  -3  -3  -1  -2  -1  -2  -1  -2  -2   2  -3 
 -1  -3  -3  -3  -1  -3  -3  -4  -3   4   2  -3   1   0  -3  -2  -1  -3  -1   3 
 -1  -2  -3  -4  -1  -2  -3  -4  -3   2   4  -2   2   0  -3  -2  -1  -2  -1   1 
 -1   2   0  -1  -3   1   1  -2  -1  -3  -2   5  -1  -3  -1   0  -1  -3  -2  -2 
 -1  -1  -2  -3  -1   0  -2  -3  -2   1   2  -1   5   0  -2  -1  -1  -1  -1   1 
 -2  -3  -3  -3  -2  -3  -3  -3  -1   0   0  -3   0   6  -4  -2  -2   1   3  -1 
 -1  -2  -2  -1  -3  -1  -1  -2  -2  -3  -3  -1  -2  -4   7  -1  -1  -4  -3  -2 
  1  -1   1   0  -1   0   0   0  -1  -2  -2   0  -1  -2  -1   4   1  -3  -2  -2 
  0  -1   0  -1  -1  -1  -1  -2  -2  -1  -1  -1  -1  -2  -1   1   5  -2  -2   0 
 -3  -3  -4  -4  -2  -2  -3  -2  -2  -3  -2  -3  -1   1  -4  -3  -2  11   2  -3 
 -2  -2  -2  -3  -2  -1  -2  -3   2  -1  -1  -2  -1   3  -3  -2  -2   2   7  -1 
  0  -3  -3  -3  -1  -2  -2  -3  -3   3   1  -2   1  -1  -2  -2   0  -3  -1   4 
  )))

(define blosum-caption "ARNDCQEGHILKMFPSTWYV")

(define (nat x) (if (= x 0) '(0) (append (nat (- x 1)) `(,x))))

(define reverse-map '())
(for-each
  (lambda (index)
    (set! reverse-map (cons (cons (string-ref blosum-caption index) index) reverse-map)))
  (nat 19))

(define (score0 so-far a b)
  (if (or (null? a) (null? b))
    so-far
    (score0
      (+ so-far (vector-ref BLOSUM62 (+ (* 20 (car a)) (car b))))
      (cdr a)
      (cdr b))))

(define magic (lambda (str) (map (lambda (char) (cdr (assoc char reverse-map))) (string->list str))))

(define (score a b) (score0 0 (magic a) (magic b)))

(define (read-fasta0 so-far port)
  (let
    ((line (read-line port)))
    (if (eof-object? line)
      so-far
      (read-fasta0 (string-append so-far line) port))))

; Read a file in Single FASTA format.
(define (read-fasta port)
  (begin
    (read-line port)
    (read-fasta0 "" port)))

(define (triples0 hash str index len)
  (if (> (+ index 3) len)
    (hash-map hash (lambda (k v) k))
    (begin
      (hash-set! hash (substring str index (+ index 3)) #t)
      (triples0
      hash
      str
      (+ index 1)
      len))))

; Harvest all triples in the given sequence.
(define (triples str) (triples0 (make-hash) str 0 (string-length str)))

(define q (call-with-input-file "titin1.fasta" read-fasta))
(define d (call-with-input-file "znf521.fasta" read-fasta))

; XXX
;(set! q (substring q 0 70))
;(set! d (substring d 0 70))

; For all triples in database,
; align the given triple (that comes from Q)
; to each of those triples in database (D)
; such that the score is above threshold.

;(define (magic10 triple d threshold so-far index)
  ;(score triple (substring d index (+ index 3))))
;(define (magic1 triple d threshold) ())

; Random alignment.
; @param a a string
; @param b a string
; @param num-trials how many random numbers are to be generated
; @return
; A list of alignments.
; An alignment is a list (i j k) that means
; a[i]..a[i+k-1] is aligned with b[j]..b[j+k-1].
(define (random-alignment0 a b max-length)
  (let*
    ((la (string-length a))
     (lb (string-length b))
     (l (min max-length la lb))
     (i (random l))
     (j (random l))
     (k (random (- l (max i j)))))
    (list i j k)))

(define (random-alignment a b max-length num-trials)
  (if (<= num-trials 0)
    '()
    (cons (random-alignment0 a b max-length) (random-alignment a b max-length (- num-trials 1)))))

; Monte Carlo Alignment.
(define (monte-carlo-alignment a b max-length min-score num-trials)
  (filter
    (lambda (alignment)
      (let*
        ((prepped (prepare-alignment a b alignment))
         (stra (car prepped))
         (strb (cdr prepped)))
        (>= (score stra strb) min-score)))
    (random-alignment a b max-length num-trials)))

; Convert integer list returned by random-alignment
; into a pair of strings that are one step closer to being displayable.
; Returns (sa . sb).  You are expected to feed this into display-alignment like this:
; (display-alignment sa sb)
(define (prepare-alignment a b l)
  (let
    ((i (first l))
     (j (second l))
     (k (third l)))
    (cons (substring a i (+ i k)) (substring b j (+ j k)))))

(define (display-alignment0 a b index len)
  (if (< index len)
    (begin
      (display (if (eq? (string-ref a index) (string-ref b index)) "|" " "))
      (display-alignment0 a b (+ index 1) len))
    (newline)))

; @param a a string.
; @param b a string.
(define (display-alignment a b i j k)
  (begin
    (display "Position A: ") (display i) (newline)
    (display "Position B: ") (display j) (newline)
    (display "Length: ") (display k) (newline)
    (display "Score: ") (display (score a b)) (newline)
    (display a) (newline)
    (display-alignment0 a b 0 (min (string-length a) (string-length b)))
    (display b) (newline)))

(for-each
  (lambda (alignment)
    (let*
      ((p (prepare-alignment q d alignment))
       (i (first alignment))
       (j (second alignment))
       (k (third alignment))
       (a (car p))
       (b (cdr p)))
      (begin
        (display-alignment a b i j k)
        (newline))))
  (monte-carlo-alignment q d 32 8 262144))

(newline)
