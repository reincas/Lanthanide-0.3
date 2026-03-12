# -*- makefile -*-
# (c) 01/2001 by Reinhard Caspary
include $(HOME)/Makefile.inc

all: wigner.ps tensors.ps hamiltonian.ps terms.ps

edit:
	xemacs Makefile *.tex &

Periodensystem.inc: Periodensystem.fig
	fig2dev -L latex -l 2 -m 0.5 Periodensystem.fig Periodensystem.inc

terms.ps: terms.tex Periodensystem.inc SE.eps

.dia.eps:
	termdia-1.3/termdia $*.dia > $*.eps

.dvi.ps:
	dvips -P600 -t a4 -t landscape -o - $*.dvi | stripland > $*.ps