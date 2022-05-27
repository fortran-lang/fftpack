# Fortran fftpack Makefile

LIB = dfftpack

FC = gfortran
FFLAGS = -O2 -fPIC

export LIB
export FC
export FFLAGS

.PHONY: build clean test

build:
	$(MAKE) -f Makefile $@ --directory=src

test: build
	$(MAKE) -f Makefile $@ --directory=test
	
bench: build
	$(MAKE) -f Makefile $@ --directory=example

clean:
	$(MAKE) -f Makefile $@ --directory=src
	$(MAKE) -f Makefile $@ --directory=test
	$(MAKE) -f Makefile $@ --directory=example