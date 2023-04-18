name=readme_regression

all: $(name).pdf

$(name).pdf: $(name).tex
	pdflatex $(name)
	pdflatex $(name)

clean:
	rm -f $~ *.aux *.idx *.log *.blg *.bbl *.toc *.dvi *.ps *.out *.gz \
			 *.fdb_latexmk *.fls *.cpt $(name).pdf
