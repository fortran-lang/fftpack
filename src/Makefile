SRCF = \
	zfftb.f90\
	cfftb1.f90\
	zfftf.f90\
	cfftf1.f90\
	zffti.f90\
	cffti1.f90\
	dcosqb.f90\
	cosqb1.f90\
	dcosqf.f90\
	cosqf1.f90\
	dcosqi.f90\
	dcost.f90\
	dcosti.f90\
	ezfft1.f90\
	dzfftb.f90\
	dzfftf.f90\
	dzffti.f90\
	passb.f90\
	passb2.f90\
	passb3.f90\
	passb4.f90\
	passb5.f90\
	passf.f90\
	passf2.f90\
	passf3.f90\
	passf4.f90\
	passf5.f90\
	radb2.f90\
	radb3.f90\
	radb4.f90\
	radb5.f90\
	radbg.f90\
	radf2.f90\
	radf3.f90\
	radf4.f90\
	radf5.f90\
	radfg.f90\
	dfftb.f90\
	rfftb1.f90\
	dfftf.f90\
	rfftf1.f90\
	dffti.f90\
	rffti1.f90\
	dsinqb.f90\
	dsinqf.f90\
	dsinqi.f90\
	dsint.f90\
	sint1.f90\
	dsinti.f90

SRCF90 = \
	fftpack.f90\
	fftpack_fft.f90\
	fftpack_ifft.f90\
	fftpack_rfft.f90\
	fftpack_irfft.f90\
	fftpack_fftshift.f90\
	fftpack_ifftshift.f90\
	fftpack_qct.f90\
	fftpack_iqct.f90\
	fftpack_dct.f90\
	rk.f90\
	fftpack_utils.f90

OBJF := $(SRCF:.f90=.o)
OBJF90 := $(SRCF90:.f90=.o)

lib$(LIB).a: $(OBJF) $(OBJF90)
	ar -rcs lib$(LIB).a $(OBJF) $(OBJF90)

shared: $(OBJ)
	$(FC) -shared -o lib$(LIB).so $(OBJ)

clean:
	rm -f -r *.o *.a *.so *.mod *.smod

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

fftpack_fft.o: fftpack.o rk.o
fftpack_ifft.o: fftpack.o rk.o
fftpack_rfft.o: fftpack.o rk.o
fftpack_irfft.o: fftpack.o rk.o
fftpack_qct.o: fftpack.o rk.o
fftpack_iqct.o: fftpack.o rk.o
fftpack_dct.o: fftpack.o rk.o
fftpack_fftshift.o: fftpack.o rk.o
fftpack_ifftshift.o: fftpack.o rk.o
fftpack_utils.o: fftpack.o

zfftb.f90: rk.o
cfftb1.f90: rk.o
zfftf.f90: rk.o
cfftf1.f90: rk.o
zffti.f90: rk.o
cffti1.f90: rk.o
dcosqb.f90: rk.o
cosqb1.f90: rk.o
dcosqf.f90: rk.o
cosqf1.f90: rk.o
dcosqi.f90: rk.o
dcost.f90: rk.o
dcosti.f90: rk.o
ezfft1.f90: rk.o
dzfftb.f90: rk.o
dzfftf.f90: rk.o
dzffti.f90: rk.o
passb.f90: rk.o
passb2.f90: rk.o
passb3.f90: rk.o
passb4.f90: rk.o
passb5.f90: rk.o
passf.f90: rk.o
passf2.f90: rk.o
passf3.f90: rk.o
passf4.f90: rk.o
passf5.f90: rk.o
radb2.f90: rk.o
radb3.f90: rk.o
radb4.f90: rk.o
radb5.f90: rk.o
radbg.f90: rk.o
radf2.f90: rk.o
radf3.f90: rk.o
radf4.f90: rk.o
radf5.f90: rk.o
radfg.f90: rk.o
dfftb.f90: rk.o
rfftb1.f90: rk.o
dfftf.f90: rk.o
rfftf1.f90: rk.o
dffti.f90: rk.o
rffti1.f90: rk.o
dsinqb.f90: rk.o
dsinqf.f90: rk.o
dsinqi.f90: rk.o
dsint.f90: rk.o
sint1.f90: rk.o
dsinti.f90: rk.o
