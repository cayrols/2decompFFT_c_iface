-include make.inc

VERBOSE=0
CP_FLAGS="-ir"
DEFAULT_PREFIX = ../
PREFIX ?= $(DEFAULT_PREFIX)

ifeq ($(VERBOSE), 1)
CP_FLAGS += "-v"
endif

.PHONY: all help install

all: help

install:
	@if [ -d $(PREFIX) ]; then \
		echo "Installation in $(PREFIX)"; \
		cp $(CP_FLAGS) examples $(PREFIX); \
		cp $(CP_FLAGS) interfaces $(PREFIX); \
		cp $(CP_FLAGS) include $(PREFIX); \
	else \
		echo "Cannot install in $(PREFIX): does not exist"; \
	fi

help:
	@echo -e "[HELP]\n"
	@echo -e "Installation of the interface by default in $(DEFAULT_PREFIX):\n\tmake install\n"
	@echo -e "To change it, either:\n"\
				"\t- make PREFIX=<2decomp_src_path> install\n"\
				"\t- create a make.inc file that contains the PREFIX variable"\
				" set to the source of the 2decomp library.\n"
	@echo "NOTE: The file make.inc must be located in the same repository as this Makefile"
