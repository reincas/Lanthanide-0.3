# -*- makefile -*-
# (c) 01/2001 by Reinhard Caspary
include $(HOME)/Makefile.inc

all: wigner.ps tensors.ps hamiltonian.ps terms.ps

edit:
	xemacs Makefile *.tex &

.dvi.ps:
	dvips -P600 -t a4 -t landscape -o - $*.dvi | stripland > $*.ps
