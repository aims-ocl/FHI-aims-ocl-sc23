
FHI-aims.pdf: *.tex FHI-aims.bbl FHI-aims.aux FHI-aims.ind
	pdflatex FHI-aims
	pdflatex FHI-aims

FHI-aims.bbl: FHI-aims.aux *.bib
	bibtex FHI-aims

FHI-aims.ind: FHI-aims.aux *.tex
	makeindex FHI-aims

FHI-aims.aux: *.tex
	pdflatex FHI-aims
	pdflatex FHI-aims


clean:
	rm -f *.aux *.bbl *.blg *.brf *.idx *.ilg *.ind *.log *.out *.toc
