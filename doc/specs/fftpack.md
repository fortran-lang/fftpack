---
title: FFTPACK
---

# The `fftpack` module

[TOC]

## `zffti`

### Description

Initializes the array `wsave` which is used in both `zfftf` and `zfftb`.   
The prime factorization of `n` together with a tabulation of the trigonometric functions are computed and
stored in `wsave`.

### Status

Experimental.

### Class

Pure function.

### Snytax

`call [[fftpack(module):zffti(interface)]](n, wsave)`

### Argument

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the sequence to be transformed.

`wsave`: Shall be a `real` array.
This argument is `intent(out)`.  
A work array which must be dimensioned at least `4*n+15`.
The same work array can be used for both `zfftf` and `zfftb`
as long as `n` remains unchanged. Different `wsave` arrays
are required for different values of `n`.   

#### Warning

The contents of `wsave` must not be changed between calls of `zfftf` or `zfftb`.

### Example

```fortran
program demo_zffti
    use fftpack, only: zffti
    complex(kind=8) :: x(4)
    real(kind=8) :: w(31)
    x = [real(kind=8) :: 1.0, 2.0, 3.0, 4.0]
    call zffti(4,w)
end program demo_zffti
```

## `zfftf`

### Description

Computes the forward complex discrete fourier
transform (the fourier analysis).   
Equivalently, `zfftf` computes
the fourier coefficients of a complex periodic sequence.
the transform is defined below at output parameter `c`.

The transform is not normalized. to obtain a normalized transform
the output must be divided by `n`. otherwise a call of `zfftf`
followed by a call of `zfftb` will multiply the sequence by `n`.

The array `wsave` which is used by subroutine `zfftf` must be
initialized by calling subroutine `zffti(n,wsave)`.

### Status

Experimental.

### Class

Pure function.

### Snytax

`call [[fftpack(module):zfftf(interface)]](n, c, wsave)`

### Argument

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the `complex` sequence `c`. The method is more efficient when `n` is the product of small primes.

`c`: Shall be a `complex` array.
This argument is `intent(inout)`.  
A `complex` array of length `n` which contains the sequence.
```
for j=1,...,n

           c(j)=the sum from k=1,...,n of

                 c(k)*exp(-i*(j-1)*(k-1)*2*pi/n)

                       where i=sqrt(-1)
```

`wsave`: Shall be a `real` array.
This argument is `intent(inout)`.  
A `real` work array which must be dimensioned at least `4n+15` in the program that calls `zfftf`. The wsave array must be initialized by calling subroutine `zffti(n,wsave)` and a different wsave array must be used for each different value of `n`. This initialization does not have to be repeated so long as `n` remains unchanged thus subsequent transforms can be obtained faster than the first. The same wsave array can be used by `zfftf` and `zfftb`.  
Contains initialization calculations which must not be destroyed between calls of subroutine `zfftf` or `zfftb`.

#### Warning

The contents of `wsave` must not be changed between calls of `zfftf` or `zfftb`.

### Example

```fortran
program demo_zfftf
    use fftpack, only: zffti, zfftf
    complex(kind=8) :: x(4)
    real(kind=8) :: w(31)
    x = [real(kind=8) :: 1.0, 2.0, 3.0, 4.0]
    call zffti(4,w)
    call zfftf(4,x,w)   !! `x` returns [(10.0,0.0), (-2.0,2.0), (-2.0,0.0), (-2.0,-2.0)].
end program demo_zfftf
```

## `zfftb`

### Description

Unnormalized inverse of `zfftf`.

Computes the backward `complex` discrete fourier
transform (the fourier synthesis). Equivalently, `zfftb` computes
a `complex` periodic sequence from its fourier coefficients.
The transform is defined below at output parameter `c`.

The transform is not normalized. to obtain a normalized transform
the output must be divided by `n`. otherwise a call of `zfftf`
followed by a call of `zfftb` will multiply the sequence by `n`.

The array `wsave` which is used by subroutine `zfftf` must be
initialized by calling subroutine `zffti(n,wsave)`.

### Status

Experimental.

### Class

Pure function.

### Snytax

`call [[fftpack(module):zfftb(interface)]](n, c, wsave)`

### Argument

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the `complex` sequence `c`. The method is more efficient when `n` is the product of small primes.

`c`: Shall be a `complex` array.
This argument is `intent(inout)`.  
A `complex` array of length `n` which contains the sequence.
```
for j=1,...,n

           c(j)=the sum from k=1,...,n of

                 c(k)*exp(-i*(j-1)*(k-1)*2*pi/n)

                       where i=sqrt(-1)
```

`wsave`: Shall be a `real` array.
This argument is `intent(inout)`.  
A `real` work array which must be dimensioned at least `4n+15` in the program that calls `zfftf`. The wsave array must be initialized by calling subroutine `zffti(n,wsave)` and a different wsave array must be used for each different value of `n`. This initialization does not have to be repeated so long as `n` remains unchanged thus subsequent transforms can be obtained faster than the first. The same wsave array can be used by `zfftf` and `zfftb`.  
Contains initialization calculations which must not be destroyed between calls of subroutine `zfftf` or `zfftb`.

#### Warning

The contents of `wsave` must not be changed between calls of `zfftf` or `zfftb`.

### Example

```fortran
program demo_zfftb
    use fftpack, only: zffti, zfftf, zfftb
    complex(kind=8) :: x(4)
    real(kind=8) :: w(31)
    x = [real(kind=8) :: 1.0, 2.0, 3.0, 4.0]
    call zffti(4,w)
    call zfftf(4,x,w)   !! `x` returns [(10.0,0.0), (-2.0,2.0), (-2.0,0.0), (-2.0,-2.0)].
    call zfftb(4,x,w)   !! `x` returns [(4.0,0.0), (8.0,0.0), (12.0,0.0), (16.0,0.0)].
end program demo_zfftb
```

## `fft`

### Description

Computes the forward complex discrete fourier
transform (the fourier analysis).   

### Status

Experimental.

### Class

Pure function.

### Snytax

`call [[fftpack(module):fft(interface)]](x [, n])`

### Argument

`x`: Shall be a `complex` array.
This argument is `intent(in)`. 

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
The needed length of the `complex` sequence `c`.

#### Warning

if `n <= size(x)`, the first `n` elements of `x` will be included in the calculation.

if `n > size(x)`, the all elements of `x` and `n-size(x)` elements filled with zeros will be included in the calculation.

### Example

```fortran
program demo_fft
    use fftpack, only: fft
    complex(kind=8) :: x(4)
    x = [real(kind=8) :: 1.0, 2.0, 3.0, 4.0]
    print *, fft(x)     !! [(10.0,0.0), (-2.0,2.0), (-2.0,0.0), (-2.0,-2.0)].
    print *, fft(x,3)   !! [(6.0,0.0), (-1.5,0.86), (-1.5,0.86)].
    print *, fft(x,5)   !! [(10.0,0.0), (-4.0,1.3), (1.5,-2.1), (1.5,2.1), (-4.0,1.3)].
end program demo_fft
```

## `ifft`

### Description

Unnormalized inverse of `fft`.

### Status

Experimental.

### Class

Pure function.

### Snytax

`call [[fftpack(module):ifft(interface)]](x [, n])`

### Argument

`x`: Shall be a `complex` array.
This argument is `intent(in)`. 

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
The needed length of the `complex` sequence `c`.

#### Warning

if `n <= size(x)`, the first `n` elements of `x` will be included in the calculation.

if `n > size(x)`, the all elements of `x` and `n-size(x)` elements filled with zeros will be included in the calculation.

### Example

```fortran
program demo_ifft
    use fftpack, only: fft, ifft
    complex(kind=8) :: x(4)
    x = [real(kind=8) :: 1.0, 2.0, 3.0, 4.0]
    print *, fft(x)             !! [(10.0,0.0), (-2.0,2.0), (-2.0,0.0), (-2.0,-2.0)].
    print *, fft(x,3)           !! [(6.0,0.0), (-1.5,0.86), (-1.5,0.86)].
    print *, fft(x,5)           !! [(10.0,0.0), (-4.0,1.3), (1.5,-2.1), (1.5,2.1), (-4.0,1.3)].
    print *, ifft(fft(x))/4.0   !! [(1.0,0.0), (2.0,0.0), (3.0,0.0), (4.0,0.0)]
    print *, ifft(fft(x), 3)    !! [(6.0,2.0), (10.3,-1.0), (13.73,-1.0)]
end program demo_ifft
```
