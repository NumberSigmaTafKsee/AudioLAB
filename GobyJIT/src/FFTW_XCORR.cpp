
std::vector<float> FFTW_CrossCorrelation(std::vector<float> & x, std::vector<float> & y)
    int M = x.size();
    int N = y.size();
    float x[M] = {0,1,0,0};
    float y[N] = {0,0,1,0,0,0,0,0};
    float in_a[M+N-1];
    std::complex<float> out_a[M+N-1];
    float in_b[M+N-1];
    std::complex<float> out_b[M+N-1];
    std::complex<float> in_rev[M+N-1];
    std::vector<float> out_rev(M+N-1);

    // Plans for forward FFTs
    fftwf_plan plan_fwd_a = fftwf_plan_dft_r2c_1d (M+N-1, in_a, reinterpret_cast<fftwf_complex*>(&out_a), FFTW_MEASURE);
    fftwf_plan plan_fwd_b = fftwf_plan_dft_r2c_1d (M+N-1, in_b, reinterpret_cast<fftwf_complex*>(&out_b), FFTW_MEASURE);

    // Plan for reverse FFT
    fftwf_plan plan_rev = fftwf_plan_dft_c2r_1d (M+N-1,reinterpret_cast<fftwf_complex*>(&in_rev), out_rev.data, FFTW_MEASURE);

    // Prepare padded input data
    std::memcpy(in_a, x.data(), sizeof(float) * M);
    std::memset(in_a + M, 0, sizeof(float) * (N-1));
    std::memset(in_b, 0, sizeof(float) * (M-1));
    std::memcpy(in_b + (M-1), y.data(), sizeof(float) * N);
     
    // Calculate the forward FFTs
    fftwf_execute(plan_fwd_a);
    fftwf_execute(plan_fwd_b);

    // Multiply in frequency domain
    for( int idx = 0; idx < M+N-1; idx++ ) {
        in_rev[idx] = std::conj(out_a[idx]) * out_b[idx]/(float)(M+N-1);
    }

    // Calculate the backward FFT
    fftwf_execute(plan_rev);

    // Clean up
    fftwf_destroy_plan(plan_fwd_a);
    fftwf_destroy_plan(plan_fwd_b);
    fftwf_destroy_plan(plan_rev);
}