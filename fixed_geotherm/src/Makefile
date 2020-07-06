.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

include Makefile.build
include Makefile.macbook

OBJECTS = \
module_definitions.o \
build_matrices.o \
closure_temps.o \
find_zc.o \
prior_zc.o \
dodson.o \
isotherms.o \
four1.o \
syn_ages3.o \
post_res.o \
posterior_dat.o \
post_full.o \
write_output.o \
initialize_parameters.o \
solve_inverse.o \
tridag_par.o \
tridag.o \
fourn.o \
kvalue.o \
clean_up.o \
glide.o 

.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90

.f.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f

.c.o:
	$(CC) $(FLAGS) $(INCLUDE) $*.c

SITE:   $(OBJECTS)
	$(F90) $(OPTIONS) -fopenmp $(OBJECTS) $(INCLUDE) $(LIBS) -o ../glide

clean:
	rm *.o *.mod

