SRCF90 = \
	bench1.f90

OBJ = $(SRCF90:%.f90=%.o) 

build: bench1.x

bench: build
	./bench1.x

bench1.x: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -L../src -l$(LIB) -I../src -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -I../src -c $<

clean:
	rm -f -r *.o *.x
