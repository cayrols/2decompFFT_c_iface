include ../../src/Makefile.inc

SRC=decomp_2d_ciface.F90

OBJ=$(SRC:.F90=.o)

all: $(OBJ)

%.o: %.F90
	@echo "Wants to compile $<"
	$(F90) $(OPTIONS) $(F90FLAGS) -I ../../include -c $<
