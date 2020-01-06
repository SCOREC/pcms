all:	pdf

dvi:	futils.dvi
pdf:	futils.pdf
ps:	futils.ps

.SUFFIXES:
.SUFFIXES: .tex .dvi .pdf .ps

.tex.dvi:
	latex $<
	@while ( grep "Rerun to get cross-references" \
	         ${<:tex=log} > /dev/null ); do \
	         latex $<; \
        done
	latex $<

.dvi.pdf:
	dvipdf $<

.dvi.ps:
	dvips $<

futils.dvi: hdfview.eps

hdfview.eps:
	convert hdfview.png hdfview.eps

clean:
	rm -f *.dvi *.aux *.log *~ *.toc

distclean:  clean
	rm -f futils.ps  hdfview.eps futils.pdf
