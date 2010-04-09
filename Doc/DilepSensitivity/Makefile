all:	dvi pdf

dvi:	dilepSensitivity.tex
	latex dilepSensitivity; latex dilepSensitivity; latex dilepSensitivity

ps:
	dvips dilepSensitivity -o dilepSensitivity.ps

pdf:
	dvipdf dilepSensitivity

clean:
	rm -rf *.aux *.dvi *.log *.pdf
