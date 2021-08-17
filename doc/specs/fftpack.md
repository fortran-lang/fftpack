---
title: FFTPACK
---

# The `fftpack` module

[TOC]

## Fourier transform of double complex periodic sequences
### `zffti`

#### Description

Initializes the array `wsave` which is used in both `zfftf` and `zfftb`.   
The prime factorization of `n` together with a tabulation of the trigonometric functions are computed and
stored in `wsave`.

#### Status

Experimental.

#### Class

Pure function.

#### Snytax

`call [[fftpack(module):zffti(interface)]](n, wsave)`

#### Argument

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the sequence to be transformed.

`wsave`: Shall be a `real` array.
This argument is `intent(out)`.  
A work array which must be dimensioned at least `4*n+15`.
The same work array can be used for both `zfftf` and `zfftb`
as long as `n` remains unchanged. Different `wsave` arrays
are required for different values of `n`.   

##### Warning

The contents of `wsave` must not be changed between calls of `zfftf` or `zfftb`.

#### Example

```fortran
program demo_zffti
    use fftpack, only: zffti
    complex(kind=8) :: x(4)
    real(kind=8) :: w(31)
    x = [real(kind=8) :: 1.0, 2.0, 3.0, 4.0]
    call zffti(4,w)
end program demo_zffti
```

### `zfftf`

#### Description

Computes the forward complex discrete fourier transform (the fourier analysis).   
Equivalently, `zfftf` computes the fourier coefficients of a complex periodic sequence.
The transform is defined below at output parameter `c`.

The transform is not normalized. To obtain a normalized transform
the output must be divided by `n`. Otherwise a call of `zfftf`
followed by a call of `zfftb` will multiply the sequence by `n`.

The array `wsave` which is used by subroutine `zfftf` must be
initialized by calling subroutine `zffti(n,wsave)`.

#### Status

Experimental.

#### Class

Pure function.

#### Snytax

`call [[fftpack(module):zfftf(interface)]](n, c, wsave)`

#### Argument

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
A `real` work array which must be dimensioned at least `4n+15` in the program that calls `zfftf`. 
The wsave array must be initialized by calling subroutine `zffti(n,wsave)` and a different `wsave` array must be used for each different value of `n`.  
This initialization does not have to be repeated so long as `n` remains unchanged thus subsequent transforms can be obtained faster than the first. 
The same `wsave` array can be used by `zfftf` and `zfftb`.  
Contains initialization calculations which must not be destroyed between calls of subroutine `zfftf` or `zfftb`.

##### Warning

The contents of `wsave` must not be changed between calls of `zfftf` or `zfftb`.

#### Example

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

### `zfftb`

#### Description

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

#### Status

Experimental.

#### Class

Pure function.

#### Snytax

`call [[fftpack(module):zfftb(interface)]](n, c, wsave)`

#### Argument

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
A `real` work array which must be dimensioned at least `4n+15` in the program that calls `zfftf`. The `wsave` array must be initialized by calling subroutine `zffti(n,wsave)` and a different `wsave` array must be used for each different value of `n`. This initialization does not have to be repeated so long as `n` remains unchanged thus subsequent transforms can be obtained faster than the first. The same `wsave` array can be used by `zfftf` and `zfftb`.  
Contains initialization calculations which must not be destroyed between calls of subroutine `zfftf` or `zfftb`.

##### Warning

The contents of `wsave` must not be changed between calls of `zfftf` or `zfftb`.

#### Example

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

### `fft`

#### Description

Computes the forward complex discrete fourier transform (the fourier analysis).   

#### Status

Experimental.

#### Class

Pure function.

#### Snytax

`result = [[fftpack(module):fft(interface)]](x [, n])`

#### Argument

`x`: Shall be a `complex` array.
This argument is `intent(in)`. 

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
The needed length of the `complex` sequence `c`.

##### Warning

if `n <= size(x)`, the first `n` elements of `x` will be included in the calculation.

if `n > size(x)`, the all elements of `x` and `n-size(x)` elements filled with zeros will be included in the calculation.

#### Return value

`fft(x)` returns the Discrete Fourier Transform (DFT) of `x` using the Fast Fourier Transform (FFT) algorithm. 

`fft(x, n)` returns `n` point DFT of `x` using the FFT algorithm.

#### Example

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

### `ifft`

#### Description

Unnormalized inverse of `fft`.

#### Status

Experimental.

#### Class

Pure function.

#### Snytax

`result = [[fftpack(module):ifft(interface)]](x [, n])`

#### Argument

`x`: Shall be a `complex` array.
This argument is `intent(in)`. 

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
The needed length of the `complex` sequence `c`.

##### Warning

if `n <= size(x)`, the first `n` elements of `x` will be included in the calculation.

if `n > size(x)`, the all elements of `x` and `n-size(x)` elements filled with zeros will be included in the calculation.

#### Return value

`ifft(x)` returns the inverse Discrete Fourier transform of `x` using the Fast Fourier Transform algorithm.. 

`ifft(x, n)` returns `n` point inverse DFT of `x` using the FFT algorithm.

#### Example

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

## Fourier transform of double real periodic sequences
### `dffti`

#### Description

Initializes the array `wsave` which is used in both `dfftf` and `dfftb`.   
The prime factorization of `n` together with a tabulation of the trigonometric functions are computed and
stored in `wsave`.

#### Status

Experimental.

#### Class

Pure function.

#### Snytax

`call [[fftpack(module):dffti(interface)]](n, wsave)`

#### Argument

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the sequence to be transformed.

`wsave`: Shall be a `real` array.
This argument is `intent(out)`.  
A work array which must be dimensioned at least `2*n+15`.
The same work array can be used for both `dfftf` and `dfftb`
as long as `n` remains unchanged. Different `wsave` arrays
are required for different values of `n`.   

##### Warning

The contents of `wsave` must not be changed between calls of `dfftf` or `dfftb`.

#### Example

```fortran
program demo_dffti
    use fftpack, only: dffti
    real(kind=8) :: x(4) = [1.0, 2.0, 3.0, 4.0]
    real(kind=8) :: w(23)
    call dffti(4,w)
end program demo_dffti
```

### `dfftf`

#### Description

Computes the fourier coefficients of a real
perodic sequence (fourier analysis). The transform is defined
below at output parameter `r`.

The transform is not normalized. To obtain a normalized transform
the output must be divided by `n`. Otherwise a call of `dfftf`
followed by a call of `dfftb` will multiply the sequence by `n`.

The array `wsave` which is used by subroutine `dfftf` must be
initialized by calling subroutine `dffti(n,wsave)`.

#### Status

Experimental.

#### Class

Pure function.

#### Snytax

`call [[fftpack(module):dfftf(interface)]](n, r, wsave)`

#### Argument

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the `real` sequence `r`. The method is more efficient when `n` is the product of small primes.
`n` may change so long as different work arrays are provided.

`r`: Shall be a `real` array.
This argument is `intent(inout)`.  
A `real` array of length `n` which contains the sequence.
```
r(1) = the sum from i=1 to i=n of r(i)

if n is even set l =n/2   , if n is odd set l = (n+1)/2

    then for k = 2,...,l

        r(2*k-2) = the sum from i = 1 to i = n of

            r(i)*cos((k-1)*(i-1)*2*pi/n)

        r(2*k-1) = the sum from i = 1 to i = n of

            -r(i)*sin((k-1)*(i-1)*2*pi/n)

if n is even

        r(n) = the sum from i = 1 to i = n of

            (-1)**(i-1)*r(i)
```

`wsave`: Shall be a `real` array.
This argument is `intent(inout)`.  
A `real` work array which must be dimensioned at least `4n+15` in the program that calls `dfftf`. 
The wsave array must be initialized by calling subroutine `dffti(n,wsave)` and a different `wsave` array must be used for each different value of `n`.  
This initialization does not have to be repeated so long as `n` remains unchanged thus subsequent transforms can be obtained faster than the first. 
The same `wsave` array can be used by `dfftf` and `dfftb`.  
Contains initialization calculations which must not be destroyed between calls of subroutine `dfftf` or `dfftb`.

##### Warning

The contents of `wsave` must not be changed between calls of `dfftf` or `dfftb`.

#### Example

```fortran
program demo_dfftf
    use fftpack, only: dffti, dfftf
    real(kind=8) :: x(4) = [1,2,3,4]
    real(kind=8) :: w(23)
    call dffti(4,w)
    call dfftf(4,x,w)   !! `x` returns [10.0, -2.0, 2.0, -2.0].
end program demo_dfftf
```

### `dfftb`

#### Description

Unnormalized inverse of `dfftf`.

Computes the backward `real` discrete fourier
transform (the fourier synthesis). Equivalently, `dfftb` computes
a `real` periodic sequence from its fourier coefficients.
The transform is defined below at output parameter `c`.

The transform is not normalized. to obtain a normalized transform
the output must be divided by `n`. otherwise a call of `dfftf`
followed by a call of `dfftb` will multiply the sequence by `n`.

The array `wsave` which is used by subroutine `dfftf` must be
initialized by calling subroutine `dffti(n,wsave)`.

#### Status

Experimental.

#### Class

Pure function.

#### Snytax

`call [[fftpack(module):dfftb(interface)]](n, r, wsave)`

#### Argument

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the `real` sequence `r`. The method is more efficient when `n` is the product of small primes.

`r`: Shall be a `real` array.
This argument is `intent(inout)`.  
A `real` array of length `n` which contains the sequence.
```
for n even and for i = 1,...,n

        r(i) = r(1)+(-1)**(i-1)*r(n)

            plus the sum from k=2 to k=n/2 of

            2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)

            -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)

for n odd and for i = 1,...,n

        r(i) = r(1) plus the sum from k=2 to k=(n+1)/2 of

            2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)

            -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
```

`wsave`: Shall be a `real` array.
This argument is `intent(inout)`.  
A `real` work array which must be dimensioned at least `2n+15` in the program that calls `dfftf`. The `wsave` array must be initialized by calling subroutine `dffti(n,wsave)` and a different `wsave` array must be used for each different value of `n`. This initialization does not have to be repeated so long as `n` remains unchanged thus subsequent transforms can be obtained faster than the first. The same `wsave` array can be used by `dfftf` and `dfftb`.  
Contains initialization calculations which must not be destroyed between calls of subroutine `dfftf` or `dfftb`.

##### Warning

The contents of `wsave` must not be changed between calls of `dfftf` or `dfftb`.

#### Example

```fortran
program demo_dfftb
    use fftpack, only: dffti, dfftf, dfftb
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    real(kind=8) :: w(31)
    call dffti(4,w)
    call dfftf(4,x,w)   !! `x` returns [10.0, -2.0, 2.0, -2.0].
    call dfftb(4,x,w)   !! `x` returns [4.0, 8.0, 12.0, 16.0], which is not normalized.
end program demo_dfftb
```

### `rfft`

#### Description

Discrete Fourier transform of a real sequence.  

#### Status

Experimental.

#### Class

Pure function.

#### Snytax

`result = [[fftpack(module):rfft(interface)]](x [, n])`

#### Argument

`x`: Shall be a `real` and rank-1 array.
This argument is `intent(in)`.  
The data to transform.

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
Defines the length of the Fourier transform. If `n` is not specified (the default) then `n = size(x)`. If `n <= size(x)`, `x` is truncated, if `n > size(x)`, `x` is zero-padded.

#### Return value

The returned real array contains:
```
[y(1),Re(y(2)),Im(y(2)),...,Re(y(n/2+1))]                if n is even
[y(1),Re(y(2)),Im(y(2)),...,Re(y(n/2+1)),Im(y(n/2+1))]   if n is odd
```
where,
```
y(j) = sum[k=1..n] x[k] * exp(-sqrt(-1)*j*k*2*pi/n)
j = 1..n
```
#### Notes

Within numerical accuracy, `y == rfft(irfft(y))/size(y)`.

#### Example

```fortran
program demo_fft
    use fftpack, only: rfft
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    print *, rfft(x,3)      !! [6.0, -1.5, 0.87].
    print *, rfft(x)        !! [10.0, -2.0, 2.0, -2.0].
    print *, rfft(x,5)      !! [10.0, -4.0, -1.3, 1.5, -2.1].
end program demo_fft
```

### `irfft`

#### Description

Unnormalized inverse of `rfft`.

#### Status

Experimental.

#### Class

Pure function.

#### Snytax

`result = [[fftpack(module):irfft(interface)]](x [, n])`

#### Argument

`x`: Shall be a `real` array.
This argument is `intent(in)`. 
Transformed data to invert.

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
Defines the length of the Fourier transform. If `n` is not specified (the default) then `n = size(x)`. If `n <= size(x)`, `x` is truncated, if `n > size(x)`, `x` is zero-padded.

#### Return value

The unnormalized inverse discrete Fourier transform.

#### Example

```fortran
program demo_irfft
    use fftpack, only: rfft, irfft
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    print *, irfft(rfft(x))/4.0   !! [1.0, 2.0, 3.0, 4.0]
    print *, irfft(rfft(x), 3)    !! [6.0, 8.53, 15.46]
end program demo_irfft
```
