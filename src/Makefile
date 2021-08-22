SRCF = \
	zfftb.f\
	cfftb1.f\
	zfftf.f\
	cfftf1.f\
	zffti.f\
	cffti1.f\
	dcosqb.f\
	cosqb1.f\
	dcosqf.f\
	cosqf1.f\
	dcosqi.f\
	dcost.f\
	dcosti.f\
	ezfft1.f\
	dzfftb.f\
	dzfftf.f\
	dzffti.f\
	passb.f\
	passb2.f\
	passb3.f\
	passb4.f\
	passb5.f\
	passf.f\
	passf2.f\
	passf3.f\
	passf4.f\
	passf5.f\
	radb2.f\
	radb3.f\
	radb4.f\
	radb5.f\
	radbg.f\
	radf2.f\
	radf3.f\
	radf4.f\
	radf5.f\
	radfg.f\
	dfftb.f\
	rfftb1.f\
	dfftf.f\
	rfftf1.f\
	dffti.f\
	rffti1.f\
	dsinqb.f\
	dsinqf.f\
	dsinqi.f\
	dsint.f\
	sint1.f\
	dsinti.f

SRCF90 = \
	fftpack.f90\
	fftpack_fft.f90\
	fftpack_ifft.f90\
	fftpack_rfft.f90\
	fftpack_irfft.f90\
	fftpack_fftshift.f90\
	fftpack_ifftshift.f90\
	fftpack_qct.f90\
	fftpack_iqct.f90

OBJF := $(SRCF:.f=.o)
OBJF90 := $(SRCF90:.f90=.o)

lib$(LIB).a: $(OBJF) $(OBJF90)
	ar -rcs lib$(LIB).a $(OBJF) $(OBJF90)

shared: $(OBJ)
	$(FC) -shared -o lib$(LIB).so $(OBJ)

clean:
	rm -f -r *.o *.a *.so *.mod *.smod

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

fftpack_fft.o: fftpack.o
fftpack_ifft.o: fftpack.o
fftpack_rfft.o: fftpack.o
fftpack_irfft.o: fftpack.o
fftpack_qct.o: fftpack.o
fftpack_iqct.o: fftpack.o
fftpack_fftshift.o: fftpack.o
fftpack_ifftshift.o: fftpack.o