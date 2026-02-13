# -*- makefile -*-
# (c) 01/2001 by Reinhard Caspary
include $(HOME)/Makefile.inc

all: test

cleanall:
	rm -f *.txt *.new

edit:
	xemacs Makefile *.py &

install:
	(cd ..; make install)

test:
	python Hss-f02.py
	python Hss-f12.py
	python Hsoo-f02.py
	python Hsoo-f12.py
	python Hci-f02.py
	python Hci-f12.py
	python Ti-f03.py

.txt.ps:
	landscape -c 160 $*.txt > $*.ps
