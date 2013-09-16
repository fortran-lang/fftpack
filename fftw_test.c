#include <time.h>
#include <fftw3.h>
int main()
{
    clock_t t1, t2, t3;
    int N;
    N = 1024*1024;
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    t1 = clock();
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    //p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_MEASURE);
    t2 = clock();

    fftw_execute(p); /* repeat as needed */
    t3 = clock();

    printf("Total time: %f ms\n", 1000. * (t3-t1) / CLOCKS_PER_SEC);
    printf("Init time : %f ms\n", 1000. * (t2-t1) / CLOCKS_PER_SEC);
    printf("Calc time : %f ms\n", 1000. * (t3-t2) / CLOCKS_PER_SEC);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
    return 0;
}
