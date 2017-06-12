all: do_simu compile_paper

do_simu:
	cd simu && make 

compile_paper:
	cd paper && make 
	cp paper/gastGaujal_AAP.pdf gastGaujal_absorbingTimeFluid_AAP2017.pdf
