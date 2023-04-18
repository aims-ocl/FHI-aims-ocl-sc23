FC := $(shell which gfortran)
FUNCS := $(sort $(wildcard *.df))


all: xc_func.f

xc_func.f: xc_func.template.f ${FUNCS} | dfauto
	cp $< $@
	FC=${FC} TMPDIR= sh $| -a -m -o $@ ${FUNCS}

dfauto:
	$(info The dfauto script can be downloaded from http://cpc.cs.qub.ac.uk/summaries/ADNY_v1_0.html.)
	$(error dfauto missing, cannot proceed)

rebuild:
	touch xc_func.template.f
	@${MAKE}

doc:
	FC=${FC} TMPDIR= sh dfauto -l -t a4paper -x 12pt *.df
	dvipdfm -o doc.pdf acg.dvi
	rm acg.*
