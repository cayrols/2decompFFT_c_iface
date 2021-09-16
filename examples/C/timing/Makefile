include ../../../src/Makefile.inc

INCLUDE = -I../../../include
LIBS = -L../../../lib -l2decomp_fft $(LIBFFT)

OBJ = 2decomp_c2c_mpi.o 
TARGETS=$(OBJ:.o=)

all: $(TARGETS)

2decomp_c2c_mpi: 2decomp_c2c_mpi.o ../../../interfaces/C/decomp_2d_ciface.o

$(TARGETS): 
	$(F90) -o $@ $^ $(LIBS)

clean:
	rm -f *.o 2decomp_c2c_mpi

%.o : %.c
	gcc $(INCLUDE) $(OPTIONS) -c $<
