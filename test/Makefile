FETCH = curl -L

SRC = \
	test_fftpack_fft.f90 \
	test_fftpack_rfft.f90 \
	test_fftpack_qct.f90 \
	test_fftpack_dct.f90 \
	test_fftpack_utils.f90 \
	test_fftpack.f90 \
	testdrive.F90
	
OBJ = $(SRC:.f90=.o)
OBJ := $(OBJ:.F90=.o)

all: tstfft \
	test_fftpack

# Orginal test
tstfft: tstfft.f
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	time ./tstfft.x

# `fftpack` fft routines
test_fftpack: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -L../src -l$(LIB) -I../src -o $@.x
	./test_fftpack.x
	
testdrive.F90:
	$(FETCH) https://github.com/fortran-lang/test-drive/raw/v0.4.0/src/testdrive.F90 > $@
	
%.o: %.F90
	$(FC) $(FFLAGS) -c $<

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
