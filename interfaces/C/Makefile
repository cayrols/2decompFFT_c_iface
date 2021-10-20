include ../../src/Makefile.inc

LIBNAME=lib2decomp_fft_iface.a
IFACE_BUILD=./build
IFACE_LIB=./lib

SRC=decomp_2d_ciface.F90

OBJ=$(addprefix $(IFACE_BUILD)/, $(SRC:.F90=.o))

.PHONY: all lib

all: $(OBJ)

lib: $(OBJ)
	@if [ ! -d $(IFACE_LIB) ]; then mkdir -pv $(IFACE_LIB); fi
	ar vr $(IFACE_LIB)/$(LIBNAME) $^
	ranlib $(IFACE_LIB)/$(LIBNAME)

$(IFACE_BUILD)/%.o: %.F90
	@echo "Wants to compile $<"
	@if [ ! -d $(IFACE_BUILD) ]; then mkdir -pv $(IFACE_BUILD); fi
	$(F90) $(OPTIONS) $(F90FLAGS) -o $@ -I ../../include -c $<

clean:
	-${RM} $(OBJ)
	-${RM} $(IFACE_LIB)/$(LIBNAME)

install:
	cp -v $(IFACE_LIB)/$(LIBNAME) ../../lib
