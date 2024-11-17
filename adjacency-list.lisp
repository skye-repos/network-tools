(load "~/quicklisp/setup.lisp")
(ql:quickload "iterate")
(use-package :iterate)

(defparameter my-list (list
						 '(2 3 4)
						 '(1 3)
						 '(2)
						 '(1 3 6)
						 '()
						 '(2 5 6)))

(defun nilsum (a b)
  (if (and (numberp a) (numberp b))
	  (+ a b) (if (numberp a)
				  a
				  b)))

(defun degree (node adj-list)
  (let* ((l1 (iter (for link in adj-list)
			   (collect (if (member node link) 1))))
		 (l2 (subst nil 0 l1)))
	(reduce #'nilsum l2)))

(defun nearest-neighbour-degree (node adj-list)
  (let* ((neighbours (nth (- node 1) adj-list)))
	(iter (for neighbour in neighbours)
	  (collect (degree neighbour adj-list)))))

(defun avg-nearest-neighbour-degree (node adj-list)
  (let* ((list (nearest-neighbour-degree node adj-list))
		 (size (length list)))
	(/ (reduce #'nilsum list) size)))
