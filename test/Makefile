FETCH = curl -L

SRCF90 = \
	test_fftpack_fft.f90 \
	test_fftpack_rfft.f90 \
	test_fftpack_utils.f90 \
	test_fftpack.f90

SRCFPP = \
	test_fftpack_qct.F90 \
	test_fftpack_dct.F90 \
	testdrive.F90

OBJ = $(SRCF90:%.f90=%.o) 
OBJ += $(SRCFPP:%.F90=%.o)

build: tstfft.x \
	test_fftpack.x
	
test: build
	./tstfft.x
	./test_fftpack.x

# Orginal test
tstfft.x: tstfft.f
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@

# `fftpack` fft routines
test_fftpack.x: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -L../src -l$(LIB) -I../src -o $@
	
testdrive.F90:
	$(FETCH) https://github.com/fortran-lang/test-drive/raw/v0.4.0/src/testdrive.F90 > $@
	
%.o: %.F90
	$(FC) $(FFLAGS) -I../src -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -I../src -c $<

test_fftpack.o: test_fftpack_fft.o \
	test_fftpack_rfft.o \
	test_fftpack_qct.o \
	test_fftpack_dct.o \
	test_fftpack_utils.o \
	testdrive.o

test_fftpack_fft.o: testdrive.o
test_fftpack_rfft.o: testdrive.o
test_fftpack_qct.o: testdrive.o
test_fftpack_dct.o: testdrive.o
test_fftpack_utils.o: testdrive.o

clean:
	rm -f -r *.o *.mod *.smod *.x testdrive.F90
