all: main.pdf

main.pdf: main.tex references.bib llncs.cls splncs04.bst
	pdflatex main
	bibtex main
	pdflatex main
	pdflatex main

clean:
	rm -f main.aux main.bbl main.blg main.log main.out main.toc main.pdf

.PHONY: all clean
