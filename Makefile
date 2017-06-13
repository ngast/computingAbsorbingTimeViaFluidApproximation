all: plot_figures compile_paper

redo_all_simu:
	cd simu && make clean && make

plot_figures:
	cd simu && make 

compile_paper:
	cd paper && make 
	cp paper/gastGaujal_AAP.pdf gastGaujal_absorbingTimeFluid_AAP2017.pdf
