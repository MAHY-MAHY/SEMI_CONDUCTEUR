# Compiler
FC = gfortran
#
#flags for debugging
FCFLAGS = -g -fbounds-check
FCFLAGS = -O2 -fopenmp
#
#flags for all
FCFLAGS += -I/usr/include
#
#libraires
LDFLAGS += -O1  -L -lblas -llapack
#
#lists of exec files
PROGRAMS = main
#
# "make" builds all
all: $(PROGRAMS)



main.o: numerics.o mod_var.o mod_tri.o allocations.o desallocations.o initialisation_syst_lin.o init_semi.o Resolution.o mod_plot.o \
	mod_new.o affiche.o   	
main  : numerics.o mod_var.o mod_tri.o allocations.o desallocations.o initialisation_syst_lin.o init_semi.o Resolution.o mod_plot.o \
	mod_new.o affiche.o #/home/joe/Documents/lapack-3.3.1/liblapack.a /home/joe/Documents/lapack-3.3.1/libblas.a






%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<




# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD		
veryclean: clean
	rm -f *~ $(PROGRAMS)

