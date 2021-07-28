# FFTPACK
A package of fortran subprograms for the fast fourier transform of periodic and other symmetric sequences.

## Getting started
### Get the code
```bash
git clone https://github.com/certik/fftpack.git
cd fftpack
```

### Build with [fortran-lang/fpm](https://github.com/fortran-lang/fpm)
Fortran Package Manager (fpm) is a great package manager and build system for Fortran.   
You can build using provided `fpm.toml`:
```bash
fpm build --flag "-O2"
fpm test --flag "-O2" tstfft
```
To use `fftpack` within your `fpm` project, add the following to your `fpm.toml` file:
```toml
[dependencies]
fftpack = { git="https://github.com/certik/fftpack.git" }
```

## Build with Make
Alternatively, you can build using provided `Makefile`:
```bash
make
```

## Links
[netlib/dfftpack1.0(fftpack4.0)](http://www.netlib.org/fftpack/)
