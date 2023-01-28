//O(N^2)
void dftProcess(const float *x, float *y, size_t N)
{
    size_t m, n;

    for (m = 0; m < N; m++) {
       y[2 * m] = 0.0f;
       y[2 * m +1] = 0.0f;

       for (n = 0; n < N; n++) {
           const float c = cosf(2 * M_PI * m * n / N);
           const float s = sinf(2 * M_PI * m * n / N);

	   y[2 * m]      += x[2 * n]      * c + x[2 * n + 1] * s;
	   y[2 * m + 1]  += x[2 * n + 1]  * c - x[2 * n] * s;
       }

    }  

}

// inverse dft
void idftProcess(const float *x, float *y, size_t N)
{
    size_t m, n;

    for (m = 0; m < N; m++) {
       y[2 * m] = 0.0f;
       y[2 * m +1] = 0.0f;

       for (n = 0; n < N; n++) {
           const float c = cosf(2 * M_PI * m * n / N);
           const float s = sinf(2 * M_PI * m * n / N);

	   y[2 * m]      += x[2 * n]      * c - x[2 * n + 1] * s;
	   y[2 * m + 1]  += x[2 * n + 1]  * c + x[2 * n] * s;
       }

       y[2 * m] = y[2 * m] / N;
       y[2 * m +1] = y[2 * m + 1] / N;
    }  

}