.SUFFIXES: .dia .eps .ps .tex .dvi

all:

clean:
	rm -f *.{bak,pyc,eps,ps,gz,log,aux,dvi}
	rm -f test.* *~

.dia.eps:
	termdia $*.dia > $*.eps

.eps.ps:
	eps2ps -s $*.eps > $*.ps

.tex.dvi:
	latex $*.tex

.dvi.ps:
	dvips -P600 -o $*.ps $*.dvi
