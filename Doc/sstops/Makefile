all:	dvi pdf

dvi:	sstops.tex
	latex sstops; latex sstops; latex sstops

ps:
	dvips sstops -o sstops.ps

pdf:
	dvipdf sstops

clean:
	rm -rf *.aux *.dvi *.log *.pdf
