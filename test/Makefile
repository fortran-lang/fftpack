all: tstfft\
fftpack_fft\
fftpack_ifft\
fftpack_rfft\
fftpack_irfft\
fftpack_fftshift\
fftpack_ifftshift

# Orginal test
tstfft: tstfft.o
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	time ./tstfft.x

# `fftpack` fft routines
fftpack_fft: test_fftpack_fft.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_fft.x

fftpack_ifft: test_fftpack_ifft.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_ifft.x

fftpack_rfft: test_fftpack_rfft.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_rfft.x

fftpack_irfft: test_fftpack_irfft.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_irfft.x

# `fftpack` utility routines
fftpack_fftshift: test_fftpack_fftshift.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_fftshift.x

fftpack_ifftshift: test_fftpack_ifftshift.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_ifftshift.x

clean:
	rm -f -r *.o *.x