# Fortran fftpack Makefile

LIB = dfftpack

FC = gfortran
FFLAGS = -O2

export LIB
export FC
export FFLAGS

.PHONY: all clean test

all:
	$(MAKE) -f Makefile --directory=src
	$(MAKE) -f Makefile --directory=test

test:
	$(MAKE) -f Makefile --directory=test

clean:
	$(MAKE) -f Makefile clean --directory=src
	$(MAKE) -f Makefile clean --directory=test