all: tstfft \
	fftpack_fft \
	fftpack_ifft \
	fftpack_rfft \
	fftpack_irfft \
	fftpack_fftshift \
	fftpack_ifftshift \
	fftpack_dzfft \
	fftpack_dcosq \
	fftpack_qct \
	fftpack_iqct \
	fftpack_dcost \
	fftpack_dct

# Orginal test
tstfft: tstfft.f
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

fftpack_dzfft: test_fftpack_dzfft.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_dzfft.x

fftpack_dcosq: test_fftpack_dcosq.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_dcosq.x

fftpack_qct: test_fftpack_qct.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_qct.x

fftpack_iqct: test_fftpack_iqct.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_iqct.x

fftpack_dcost: test_fftpack_dcost.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_dcost.x

fftpack_dct: test_fftpack_dct.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_dct.x

# `fftpack` utility routines
fftpack_fftshift: test_fftpack_fftshift.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_fftshift.x

fftpack_ifftshift: test_fftpack_ifftshift.f90
	$(FC) $(FFLAGS) $< -L../src -l$(LIB) -I../src -o $@.x
	./fftpack_ifftshift.x

clean:
	rm -f -r *.o *.x
