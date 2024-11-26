.PHONY: default

default:
	make lisp
	make man.pdf

lisp: lisp.c
	gcc -Wall -O0 -o lisp lisp.c

man.pdf: man.tex
	pdflatex man.tex
