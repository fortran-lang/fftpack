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

Pure subroutine.

#### Syntax

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
    complex(kind=8) :: x(4) = [1.0, 2.0, 3.0, 4.0]
    real(kind=8) :: w(31)
    call zffti(4,w)
end program demo_zffti
```

### `zfftf`

#### Description

Computes the forward complex discrete fourier transform (the fourier analysis).   
Equivalently, `zfftf` computes the fourier coefficients of a complex periodic sequence.
The transform is defined below at output parameter `c`.

The transform is not normalized. To obtain a normalized transform the output must be divided by `n`. Otherwise a call of `zfftf` followed by a call of `zfftb` will multiply the sequence by `n`.

The array `wsave` which is used by subroutine `zfftf` must be initialized by calling subroutine `zffti(n,wsave)`.

#### Status

Experimental.

#### Class

Pure subroutine.

#### Syntax

`call [[fftpack(module):zfftf(interface)]](n, c, wsave)`

#### Argument

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the `complex` sequence `c`. The method is more efficient when `n` is the product of small primes.

`c`: Shall be a `complex` and rank-1 array.
This argument is `intent(inout)`.  
A `complex` array of length `n` which contains the sequence.
```
for j=1,...,n

           c(j)=the sum from k=1,...,n of

                 c(k)*exp(-i*(j-1)*(k-1)*2*pi/n)

                       where i=sqrt(-1)
```

`wsave`: Shall be a `real` array.
This argument is `intent(in)`.  
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

Computes the backward `complex` discrete fourier transform (the fourier synthesis). 
Equivalently, `zfftb` computes a `complex` periodic sequence from its fourier coefficients. 
The transform is defined below at output parameter `c`.

The transform is not normalized. to obtain a normalized transform the output must be divided by `n`. Otherwise a call of `zfftf` followed by a call of `zfftb` will multiply the sequence by `n`.

The array `wsave` which is used by subroutine `zfftf` must be initialized by calling subroutine `zffti(n,wsave)`.

#### Status

Experimental.

#### Class

Pure subroutine.

#### Syntax

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
This argument is `intent(in)`.  
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

#### Syntax

`result = [[fftpack(module):fft(interface)]](x [, n])`

#### Argument

`x`: Shall be a `complex` and rank-1 array.
This argument is `intent(in)`. 

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
Defines the length of the Fourier transform. If `n` is not specified (the default) then `n = size(x)`. If `n <= size(x)`, `x` is truncated, if `n > size(x)`, `x` is zero-padded.

#### Return value

Returns a `complex` and rank-1 array, the Discrete Fourier Transform (DFT) of `x`.

#### Notes

Within numerical accuracy, `x == ifft(fft(x))/size(x)`.

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

#### Syntax

`result = [[fftpack(module):ifft(interface)]](x [, n])`

#### Argument

`x`: Shall be a `complex` and rank-1 array.
This argument is `intent(in)`. 

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
Defines the length of the Fourier transform. If `n` is not specified (the default) then `n = size(x)`. If `n <= size(x)`, `x` is truncated, if `n > size(x)`, `x` is zero-padded.

#### Return value

Returns a `complex` and rank-1 array, the unnormalized inverse Discrete Fourier Transform (DFT) of `x`.

#### Example

```fortran
program demo_ifft
    use fftpack, only: fft, ifft
    complex(kind=8) :: x(4) = [1.0, 2.0, 3.0, 4.0]
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

Pure subroutine.

#### Syntax

`call [[fftpack(module):dffti(interface)]](n, wsave)`

#### Argument

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the sequence to be transformed.

`wsave`: Shall be a `real` array.
This argument is `intent(out)`.  
A work array which must be dimensioned at least `2*n+15`.
The same work array can be used for both `dfftf` and `dfftb` as long as `n` remains unchanged. 
Different `wsave` arrays are required for different values of `n`.   

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

Computes the fourier coefficients of a real perodic sequence (fourier analysis). 
The transform is defined below at output parameter `r`.

The transform is not normalized. To obtain a normalized transform the output must be divided by `n`. Otherwise a call of `dfftf` followed by a call of `dfftb` will multiply the sequence by `n`.

The array `wsave` which is used by subroutine `dfftf` must be initialized by calling subroutine `dffti(n,wsave)`.

#### Status

Experimental.

#### Class

Pure subroutine.

#### Syntax

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
This argument is `intent(in)`.  
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
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    real(kind=8) :: w(23)
    call dffti(4,w)
    call dfftf(4,x,w)   !! `x` returns [10.0, -2.0, 2.0, -2.0].
end program demo_dfftf
```

### `dfftb`

#### Description

Unnormalized inverse of `dfftf`.

Computes the backward `real` discrete fourier transform (the fourier synthesis). 
Equivalently, `dfftb` computes a `real` periodic sequence from its fourier coefficients.
The transform is defined below at output parameter `c`.

The transform is not normalized. To obtain a normalized transform the output must be divided by `n`. Otherwise a call of `dfftf` followed by a call of `dfftb` will multiply the sequence by `n`.

The array `wsave` which is used by subroutine `dfftf` must be initialized by calling subroutine `dffti(n,wsave)`.

#### Status

Experimental.

#### Class

Pure subroutine.

#### Syntax

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
This argument is `intent(in)`.  
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

#### Syntax

`result = [[fftpack(module):rfft(interface)]](x [, n])`

#### Argument

`x`: Shall be a `real` and rank-1 array.
This argument is `intent(in)`.  
The data to transform.

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
Defines the length of the Fourier transform. If `n` is not specified (the default) then `n = size(x)`. If `n <= size(x)`, `x` is truncated, if `n > size(x)`, `x` is zero-padded.

#### Return value

Returns a `real` and rank-1 array, the Discrete Fourier Transform (DFT) of `x`.

#### Notes

Within numerical accuracy, `y == rfft(irfft(y))/size(y)`.

#### Example

```fortran
program demo_rfft
    use fftpack, only: rfft
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    print *, rfft(x,3)      !! [6.0, -1.5, 0.87].
    print *, rfft(x)        !! [10.0, -2.0, 2.0, -2.0].
    print *, rfft(x,5)      !! [10.0, -4.0, -1.3, 1.5, -2.1].
end program demo_rfft
```

### `irfft`

#### Description

Unnormalized inverse of `rfft`.

#### Status

Experimental.

#### Class

Pure function.

#### Syntax

`result = [[fftpack(module):irfft(interface)]](x [, n])`

#### Argument

`x`: Shall be a `real` array.
This argument is `intent(in)`. 
Transformed data to invert.

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
Defines the length of the Fourier transform. If `n` is not specified (the default) then `n = size(x)`. If `n <= size(x)`, `x` is truncated, if `n > size(x)`, `x` is zero-padded.

#### Return value

Returns a `real` and rank-1 array, the unnormalized inverse discrete Fourier transform.

#### Example

```fortran
program demo_irfft
    use fftpack, only: rfft, irfft
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    print *, irfft(rfft(x))/4.0   !! [1.0, 2.0, 3.0, 4.0]
    print *, irfft(rfft(x), 3)    !! [6.0, 8.53, 15.46]
end program demo_irfft
```

## Simplified fourier transform of double real periodic sequences

### `dzffti`

#### Description

Initializes the array `wsave` which is used in both `dzfftf` and `dzfftb`. 
The prime factorization of `n` together with a tabulation of the trigonometric functions are computed and stored in `wsave`.

#### Status

Experimental

#### Class

Prue function.

#### Syntax

`call [[fftpack(module):dzffti(interface)]](n, wsave)`

#### Arguments

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the sequence to be transformed.

`wsave`: Shall be a `real` and rank-1 array.
This argument is `intent(out)`.  
A work array which must be dimensioned at least `3*n+15`.
The same work array can be used for both `dzfftf` and `dzfftb` as long as n remains unchanged. 
Different `wsave` arrays are required for different values of `n`.

##### Warning

The contents of `wsave` must not be changed between calls of `dzfftf` or `dzfftb`.

#### Example

```fortran
program demo_dzffti
    use fftpack, only: dzffti
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    real(kind=8) :: w(3*4 + 15)
    call dzffti(4, w)   !! Initializes the array `w` which is used in both `dzfftf` and `dzfftb`. 
end program demo_dzffti
```

### `dzfftf`

#### Description

Computes the fourier coefficients of a `real` perodic sequence (fourier analysis). 
The transform is defined below at output parameters `azero`, `a` and `b`. 
`dzfftf` is a simplified but **slower version** of `dfftf`.

#### Status

Experimental

#### Class

Pure subroutine.

#### Syntax

`call [[fftpack(module):dzfftf(interface)]](n, r, azero, a, b, wsave)`

#### Arguments

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the array `r` to be transformed.  
The method is most efficient when `n` is the product of small primes.

`r`: Shall be a `real` and rank-1 array.
This argument is `intent(in)`.  
A `real` array of length `n` which contains the sequence to be transformed. `r` is not destroyed.

`azero`: Shall be a `real` scalar.
This argument is `intent(out)`.  
The sum from `i=1` to `i=n` of `r(i)/n`.

`a`, `b`: Shall be a `real` and rank-1 array.
This argument is `intent(out)`.  
```
for n even b(n/2)=0. and a(n/2) is the sum from i=1 to i=n of (-1)**(i-1)*r(i)/n

for n even define kmax=n/2-1
for n odd  define kmax=(n-1)/2

then for  k=1,...,kmax

        a(k) equals the sum from i=1 to i=n of

            2./n*r(i)*cos(k*(i-1)*2*pi/n)

        b(k) equals the sum from i=1 to i=n of

            2./n*r(i)*sin(k*(i-1)*2*pi/n)
```

`wsave`: Shall be a `real` and rank-1 array.
This argument is `intent(in)`.
A work array which must be dimensioned at least `3*n+15`. 
In the program that calls `dzfftf`. The `wsave` array must be initialized by calling subroutine `dzffti(n,wsave)` and a different `wsave` array must be used for each different value of `n`. 
This initialization does not have to be repeated so long as `n` remains unchanged thus subsequent transforms can be obtained faster than the first. 
The same `wsave` array can be used by `dzfftf` and `dzfftb`.

#### Example

```fortran
program demo_dzfftf
    use fftpack, only: dzffti, dzfftf
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    real(kind=8) :: w(3*4 + 15)
    real(kind=8) :: azero, a(4/2), b(4/2)
    call dzffti(4, w)
    call dzfftf(4, x, azero, a, b, w)   !! `azero`: 2.5; `a`: [-1.0, -0.5]; `b`: [-1.0, -0.0]
end program demo_dzfftf
```

### `dzfftb`

#### Description

Computes a `real` perodic sequence from its fourier coefficients (fourier synthesis). 
The transform is defined below at output parameter `r`. 
`dzfftb` is a simplified but **slower version** of `dfftb`.

#### Status

Experimental

#### Class

Pure subroutine.

#### Syntax

`call [[fftpack(module):dzfftb(interface)]](n, r, azero, a, b, wsave)`

#### Arguments

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the output array `r`.  
The method is most efficient when `n` is the product of small primes.

`r`: Shall be a `real` and rank-1 array.
This argument is `intent(out)`.  
```
if n is even define kmax=n/2
if n is odd  define kmax=(n-1)/2

then for i=1,...,n

        r(i)=azero plus the sum from k=1 to k=kmax of

        a(k)*cos(k*(i-1)*2*pi/n)+b(k)*sin(k*(i-1)*2*pi/n)
```
Complex notation:
```
for j=1,...,n

r(j) equals the sum from k=-kmax to k=kmax of

        c(k)*exp(i*k*(j-1)*2*pi/n)

where

        c(k) = .5*cmplx(a(k),-b(k))   for k=1,...,kmax

        c(-k) = conjg(c(k))

        c(0) = azero

            and i=sqrt(-1)
```
Amplitude - phase notation:
```
for i=1,...,n

r(i) equals azero plus the sum from k=1 to k=kmax of

        alpha(k)*cos(k*(i-1)*2*pi/n+beta(k))

where

        alpha(k) = sqrt(a(k)*a(k)+b(k)*b(k))

        cos(beta(k))=a(k)/alpha(k)

        sin(beta(k))=-b(k)/alpha(k)
```

`azero`: Shall be a `real` scalar.
This argument is `intent(in)`.  
The constant fourier coefficient.

`a`, `b`: Shall be a `real` and rank-1 array.
This argument is `intent(in)`.  
Arrays which contain the remaining fourier coefficients these arrays are not destroyed.
The length of these arrays depends on whether `n` is even or odd.
```
if n is even n/2    locations are required
if n is odd (n-1)/2 locations are required
```

`wsave`: Shall be a `real` and rank-1 array.
This argument is `intent(in)`.
A work array which must be dimensioned at least `3*n+15`. 
In the program that calls `dzfftf`. The `wsave` array must be initialized by calling subroutine `dzffti(n,wsave)` and a different `wsave` array must be used for each different value of `n`. 
This initialization does not have to be repeated so long as `n` remains unchanged thus subsequent transforms can be obtained faster than the first. 
The same `wsave` array can be used by `dzfftf` and `dzfftb`.

#### Example

```fortran
program demo_dzfftb
    use fftpack, only: dzffti, dzfftf, dzfftb
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    real(kind=8) :: w(3*4 + 15)
    real(kind=8) :: azero, a(4/2), b(4/2)
    call dzffti(4, w)
    call dzfftf(4, x, azero, a, b, w)   !! `azero`: 2.5; `a`: [-1.0, -0.5]; `b`: [-1.0, -0.0]
    x = 0.0
    call dzfftb(4, x, azero, a, b, w)   !! `x`: [1.0, 2.0, 3.0, 4.0]
end program demo_dzfftb
```

## Cosine transform with odd wave numbers

### `dcosqi`

#### Description

Initializes the array `wsave` which is used in both `dcosqf` and `dcosqb`. 
The prime factorization of `n` together with
a tabulation of the trigonometric functions are computed and
stored in `wsave`.

#### Status

Experimental

#### Class

Pure subroutine.

#### Syntax

`call [[fftpack(module):dcosqi(interface)]](n, wsave)`

#### Arguments

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the array to be transformed. 
The method is most efficient when `n` is a product of small primes.

`wsave`: Shall be a `real` and rank-1 array.
This argument is `intent(out)`.  
A work array which must be dimensioned at least `3*n+15`.
The same work array can be used for both `dcosqf` and `dcosqb`
as long as `n` remains unchanged. 
Different `wsave` arrays are required for different values of `n`.
The contents of `wsave` must not be changed between calls of `dcosqf` or `dcosqb`.

#### Example

```fortran
program demo_dcosqi
    use fftpack, only: dcosqi
    real(kind=8) :: w(3*4 + 15)
    call dcosqi(4, w)   !! Initializes the array `w` which is used in both `dcosqf` and `dcosqb`. 
end program demo_dcosqi
```

### `dcosqf`

#### Decsription

Computes the fast fourier transform of quarter wave data. 
That is, `dcosqf` computes the coefficients in a cosine series representation with only odd wave numbers. 
The transform is defined below at output parameter `x`.

`dcosqf` is the unnormalized inverse of `dcosqb` since a call of `dcosqf` followed by a call of `dcosqb` will multiply the input sequence `x` by `4*n`.

The array `wsave` which is used by subroutine `dcosqf` must be initialized by calling subroutine `dcosqi(n,wsave)`.

#### Status

Experimental

#### Class

Pure subroutine.

#### Syntax

`call [[fftpack(module):dcosqf(interface)]](n, x, wsave)`

#### Arguments

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the array `x` to be transformed. 
The method is most efficient when `n` is a product of small primes.

`x`: Shall be a `real` and rank-1 array.
This argument is `intent(inout)`.  
An array which contains the sequence to be transformed.
```
for i=1,...,n

        x(i) = x(1) plus the sum from k=2 to k=n of

        2*x(k)*cos((2*i-1)*(k-1)*pi/(2*n))

        a call of dcosqf followed by a call of
        cosqb will multiply the sequence x by 4*n.
        therefore dcosqb is the unnormalized inverse
        of dcosqf.
```

`wsave`: Shall be a `real` and rank-1 array.
This argument is `intent(in)`.  
A work array which must be dimensioned at least `3*n+15`
in the program that calls `dcosqf`. 
The `wsave` array must be initialized by calling subroutine `dcosqi(n,wsave)` and a different `wsave` array must be used for each different value of `n`. 
This initialization does not have to be repeated so long as `n` remains unchanged thus subsequent transforms can be obtained faster than the first.

##### Warning

`wsave` contains initialization calculations which must not be destroyed between calls of `dcosqf` or `dcosqb`.

#### Example

```fortran
program demo_dcosqf
    use fftpack, only: dcosqi, dcosqf
    real(kind=8) :: w(3*4 + 15)
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    call dcosqi(4, w) 
    call dcosqf(4, x, w)    !! `x`: [12.0, -9.10, 2.62, -1.51]
end program demo_dcosqf
```

### `dcosqb`

#### Decsription

Computes the fast fourier transform of quarter wave data.
That is, `dcosqb` computes a sequence from its representation in terms of a cosine series with odd wave numbers. 
The transform is defined below at output parameter `x`.

`dcosqb` is the unnormalized inverse of `dcosqf` since a call of `dcosqb` followed by a call of `dcosqf` will multiply the input sequence `x` by `4*n`.

The array `wsave` which is used by subroutine `dcosqb` must be initialized by calling subroutine `dcosqi(n,wsave)`.

#### Status

Experimental

#### Class

Pure subroutine.

#### Syntax

`call [[fftpack(module):dcosqf(interface)]](n, x, wsave)`

#### Arguments

`n`: Shall be an `integer` scalar.
This argument is `intent(in)`.  
The length of the array `x` to be transformed. 
The method is most efficient when `n` is a product of small primes.

`x`: Shall be a `real` and rank-1 array.
This argument is `intent(inout)`.  
An array which contains the sequence to be transformed.
```
for i=1,...,n

        x(i)= the sum from k=1 to k=n of

        4*x(k)*cos((2*k-1)*(i-1)*pi/(2*n))

        a call of dcosqb followed by a call of
        dcosqf will multiply the sequence x by 4*n.
        therefore dcosqf is the unnormalized inverse
        of dcosqb.
```

`wsave`: Shall be a `real` and rank-1 array.
This argument is `intent(in)`.  
A work array which must be dimensioned at least `3*n+15`
in the program that calls `dcosqb`. 
The `wsave` array must be initialized by calling subroutine `dcosqi(n,wsave)` and a different `wsave` array must be used for each different value of `n`. 
This initialization does not have to be repeated so long as `n` remains unchanged thus subsequent transforms can be obtained faster than the first.

##### Warning

`wsave` contains initialization calculations which must not be destroyed between calls of `dcosqf` or `dcosqb`.

#### Example

```fortran
program demo_dcosqb
    use fftpack, only: dcosqi, dcosqf, dcosqb
    real(kind=8) :: w(3*4 + 15)
    real(kind=8) :: x(4) = [4, 3, 5, 10]
    call dcosqi(4, w) 
    call dcosqf(4, x, w) 
    call dcosqb(4, x, w)    !! `x`: [1.0, 2.0, 3.0, 4.0] * 4 * n, n = 4, which is unnormalized.
end program demo_dcosqb
```

### `qct`

#### Description

Forward transform of quarter wave data.

#### Status

Experimental.

#### Class

Pure function.

#### Syntax

`result = [[fftpack(module):qct(interface)]](x [, n])`

#### Argument

`x`: Shall be a `real` and rank-1 array.
This argument is `intent(in)`.  
The data to transform.

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
Defines the length of the Fourier transform. If `n` is not specified (the default) then `n = size(x)`. If `n <= size(x)`, `x` is truncated, if `n > size(x)`, `x` is zero-padded.

#### Return value

Returns a `real` and rank-1 array, the Quarter-Cosine Transform (QCT) of `x`.

#### Notes

Within numerical accuracy, `x == iqct(qct(x))/(4*size(x))`.

#### Example

```fortran
program demo_qct
    use fftpack, only: qct
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    print *, qct(x,3)      !! [7.4, -5.0, 0.53].
    print *, qct(x)        !! [12.0, -9.10, 2.62, -1.51].
    print *, qct(x,5)      !! [14.4, -6.11, -5.0, 4.4, -2.65].
end program demo_qct
```

### `iqct`

#### Description

Unnormalized inverse of `qct`.

#### Status

Experimental.

#### Class

Pure function.

#### Syntax

`result = [[fftpack(module):iqct(interface)]](x [, n])`

#### Argument

`x`: Shall be a `real` array.
This argument is `intent(in)`. 
Transformed data to invert.

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
Defines the length of the Fourier transform. If `n` is not specified (the default) then `n = size(x)`. If `n <= size(x)`, `x` is truncated, if `n > size(x)`, `x` is zero-padded.

#### Return value

Returns a `real` and rank-1 array, the unnormalized inverse Quarter-Cosine Transform.

#### Example

```fortran
program demo_iqct
    use fftpack, only: qct, iqct
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    print *, iqct(qct(x))/(4.0*4.0)         !! [1.0, 2.0, 3.0, 4.0]
    print *, iqct(qct(x), 3)/(4.0*3.0)      !! [1.84, 2.71, 5.47]
end program demo_iqct
```

## Cosine transform of a real even sequence

### `dcosti`

#### Description

Initializes the array `wsave` which is used in subroutine `dcost`. 
The prime factorization of `n` together with a tabulation of the trigonometric functions are computed and stored in `wsave`.

#### Status

Experimental

#### Class

Pure subroutine.

#### Syntax

`call [[fftpack(module):dcosti(interface)]](n , wsave)`

#### Arguments

`n`: Shall be a `integer` scalar.
This argument is `intent(in)`.  
The length of the sequence to be transformed.  
The method is most efficient when n-1 is a product of small primes.

`wsave`: Shall be a `real` and rank-1 array.
This argument is `intent(out)`.  
A work array which must be dimensioned at least `3*n+15`. 
Different `wsave` arrays are required for different values of `n`. 
The contents of `wsave` must not be changed between calls of `dcost`.

#### Example

```fortran
program demo_dcosti
    use fftpack, only: dcosti
    real(kind=8) :: w(3*4 + 15)
    call dcosti(4, w)   !! Initializes the array `w` which is used in subroutine `dcost`. 
end program demo_dcosti
```

### `dcost`

#### Description

Computes the discrete fourier cosine transform of an even sequence `x(i)`. 
The transform is defined below at output parameter `x`.

`dcost` is the unnormalized inverse of itself since a call of `dcost` followed by another call of `dcost` will multiply the input sequence `x` by `2*(n-1)`. 
The transform is defined below at output parameter `x`.

The array `wsave` which is used by subroutine `dcost` must be initialized by calling subroutine `dcosti(n,wsave)`.

#### Status

Experimental

#### Class

Pure subroutine.

#### Syntax

`call [[fftpack(module):dcost(interface)]](n, x, wsave)`

#### Arguments

`n`: Shall be a `integer` scalar.
This argument is `intent(in)`.  
The length of the sequence `x`. 
`n` must be greater than `1`.
The method is most efficient when `n-1` is a product of small primes.

`x`: Shall be a `real` and rank-1 array.
This argument is `intent(inout)`.
An array which contains the sequence to be transformed.
```
for i=1,...,n

    x(i) = x(1)+(-1)**(i-1)*x(n)

        + the sum from k=2 to k=n-1

            2*x(k)*cos((k-1)*(i-1)*pi/(n-1))

        a call of dcost followed by another call of
        dcost will multiply the sequence x by 2*(n-1)
        hence dcost is the unnormalized inverse
        of itself.
```

`wsave`: Shall be a `real` and rank-1 array.
This argument is `intent(in)`.  
A work array which must be dimensioned at least `3*n+15` in the program that calls `dcost`. 
The `wsave` array must be initialized by calling subroutine `dcosti(n,wsave)` and a different `wsave` array must be used for each different value of `n`. 
This initialization does not have to be repeated so long as `n` remains unchanged thus subsequent
transforms can be obtained faster than the first.
Contains initialization calculations which must not be destroyed between calls of `dcost`.

#### Example

```fortran
program demo_dcost
    use fftpack, only: dcosti, dcost
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    real(kind=8) :: w(3*4 + 15)
    call dcosti(4, w)
    call dcost(4, x, w)     !! Computes the discrete fourier cosine (forward) transform of an even sequence, `x`(unnormalized): [15.0, -4.0, 0.0, -1.0]
    call dcost(4, x, w)     !! Computes the discrete fourier cosine (backward) transform of an even sequence, `x`(unnormalized): [6.0, 12.0, 18.0, 24.0]
end program demo_dcost
```

### `dct`

#### Description

Discrete fourier cosine (forward) transform of an even sequence. 

#### Status

Experimental.

#### Class

Pure function.

#### Syntax

`result = [[fftpack(module):dct(interface)]](x [, n])`

#### Argument

`x`: Shall be a `real` and rank-1 array.
This argument is `intent(in)`.  
The data to transform.

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
Defines the length of the Fourier transform. If `n` is not specified (the default) then `n = size(x)`. If `n <= size(x)`, `x` is truncated, if `n > size(x)`, `x` is zero-padded.

#### Return value

Returns a `real` and rank-1 array, the Discrete-Cosine Transform (DCT) of `x`.

#### Notes

Within numerical accuracy, `y == dct(idct(y))/2*(size(y) - 1)`.

#### Example

```fortran
program demo_dct
    use fftpack, only: dct
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    print *, dct(x,3)      !! [8.0, -2.0, 0.0].
    print *, dct(x)        !! [15.0, -4.0, 0.0, -1.0].
    print *, dct(x,5)      !! [19.0, -1.8, -5.0, 3.8, -5.0].
    print *, dct(dct(x))/(2*(4 - 1))   !! (normalized): [1.0, 2.0, 3.0, 4.0]
end program demo_dct
```

### `idct`

#### Description

Unnormalized inverse of `dct`.
In fact, `idct` and `dct` have the same effect, `idct` = `dct`.

#### Status

Experimental.

#### Class

Pure function.

#### Syntax

`result = [[fftpack(module):idct(interface)]](x [, n])`

#### Argument

`x`: Shall be a `real` array.
This argument is `intent(in)`. 
Transformed data to invert.

`n`: Shall be an `integer` scalar.
This argument is `intent(in)` and `optional`.  
Defines the length of the Fourier transform. If `n` is not specified (the default) then `n = size(x)`. If `n <= size(x)`, `x` is truncated, if `n > size(x)`, `x` is zero-padded.

#### Return value

Returns a `real` and rank-1 array, the inverse Discrete-Cosine Transform (iDCT) of `x`.

#### Example

```fortran
program demo_idct
    use fftpack, only: dct, idct
    real(kind=8) :: x(4) = [1, 2, 3, 4]
    print *, idct(dct(x))/(2*(4-1))   !! (normalized):   [1.0, 2.0, 3.0, 4.0]
    print *, idct(idct(x))/(2*(4-1))  !! (normalized):   [1.0, 2.0, 3.0, 4.0]
    print *, idct(dct(x), 3)          !! (unnormalized): [7.0, 15.0, 23.0]
end program demo_idct
```


## Utility functions

### `fftshift`

#### Description

Rearranges the Fourier transform by moving the zero-frequency component to the center of the array.  

#### Status

Experimental.

#### Class

Pure function.

#### Syntax

`result = [[fftpack(module):fftshift(interface)]](x)`

#### Argument

`x`: Shall be a `complex/real` and rank-1 array.
This argument is `intent(in)`. 

#### Return value

Returns the `complex/real` and rank-1 Fourier transform by moving the zero-frequency component to the center of the array.

#### Example

```fortran
program demo_fftshift
    use fftpack, only: fftshift
    complex(kind=8) :: c(5) = [1, 2, 3, 4, 5]
    real(kind=8) :: x(5) = [1, 2, 3, 4, 5]
    print *, fftshift(c(1:4))   !! [(3.0,0.0), (4.0,0.0), (1.0,0.0), (2.0,0.0)]
    print *, fftshift(c)        !! [(4.0,0.0), (5.0,0.0), (1.0,0.0), (2.0,0.0), (3.0,0.0)]
    print *, fftshift(x(1:4))   !! [3.0, 4.0, 1.0, 2.0]
    print *, fftshift(x)        !! [4.0, 5.0, 1.0, 2.0, 3.0]
end program demo_fftshift
```

### `ifftshift`

#### Description

Rearranges the Fourier transform with zero frequency shifting back to the original transform output. In other words, `ifftshift` is the result of undoing `fftshift`.

#### Status

Experimental.

#### Class

Pure function.

#### Syntax

`result = [[fftpack(module):ifftshift(interface)]](x)`

#### Argument

`x`: Shall be a `complex/real` and rank-1 array.
This argument is `intent(in)`. 

#### Return value

Returns the `complex/real` and rank-1 Fourier transform with zero frequency shifting back to the original transform output.

#### Example

```fortran
program demo_ifftshift
    use fftpack, only: fftshift, ifftshift
    complex(kind=8) :: c(5) = [1, 2, 3, 4, 5]
    real(kind=8) :: x(5) = [1, 2, 3, 4, 5]
    print *, ifftshift(fftshift(c(1:4)))   !! [(1.0,0.0), (2.0,0.0), (3.0,0.0), (4.0,0.0)]
    print *, ifftshift(fftshift(c) )       !! [(1.0,0.0), (2.0,0.0), (3.0,0.0), (4.0,0.0), (5.0,0.0)]
    print *, ifftshift(fftshift(x(1:4)))   !! [1.0, 2.0, 3.0, 4.0]
    print *, ifftshift(fftshift(x))        !! [1.0, 2.0, 3.0, 4.0, 5.0]
end program demo_ifftshift
```
