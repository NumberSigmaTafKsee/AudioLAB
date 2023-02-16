#pragma once

namespace AudioDSP
{

    struct FFTCorrelationDouble
    {
        FFTPlanRealDouble fftPlan1,fftPlan2;
        size_t  fftSize;
         
        std::vector<double> t1,t2;
        std::vector<std::complex<double>> o1,o2;
        
        FFTCorrelationDouble(size_t blocks)
        {           
            fftSize = blocks;
            
            fftPlan1.init(fftSize);            
            fftPlan2.init(fftSize);            
            
            t1.resize(fftSize);
            t2.resize(fftSize);
            o1.resize(fftSize);
            o2.resize(fftSize);
            
        }

        void ProcessBlock(size_t n, const double * in1, const double * in2, std::complex<double> * out)
        {
            memset(t1.data(),0,t1.size()*sizeof(double));
            memset(t2.data(),0,t2.size()*sizeof(double));

            memcpy(t1.data(),in1,n*sizeof(double));
            memcpy(t2.data(),in2,n*sizeof(double));
            
            fft(fftPlan1,t1.data(),o1.data(),false);
            fft(fftPlan1,t2.data(),o2.data(),false);
                                    
            for(size_t i = 1; i < t1.size()/2-1; i++)
            {
                out[i] = std::conj(o1[i]) * o2[i];
            }                        
        }
    };


    struct FFTCorrelationFloat
    {
        FFTPlanRealFloat fftPlan1,fftPlan2;
        size_t  fftSize;
         
        std::vector<float> t1,t2;
        std::vector<std::complex<float>> o1,o2;
        
        FFTCorrelationFloat(size_t blocks)
        {           
            fftSize = blocks;
            fftPlan1.init(fftSize);            
            fftPlan2.init(fftSize);            
            
            t1.resize(fftSize);
            t2.resize(fftSize);
            o1.resize(fftSize);
            o2.resize(fftSize);
        }

        void ProcessBlock(size_t n, const float* in1, const float* in2, std::complex<float> * out)
        {
            memset(t1.data(),0,t1.size()*sizeof(float));
            memset(t2.data(),0,t2.size()*sizeof(float));
            
            memcpy(t1.data(),in1,n*sizeof(float));
            memcpy(t2.data(),in2,n*sizeof(float));
            
            fft(fftPlan1,t1.data(),o1.data(),false);
            fft(fftPlan1,t2.data(),o2.data(),false);
                                    
            for(size_t i = 1; i < t1.size()/2-1; i++)
            {
                out[i] = std::conj(o1[i]) * o2[i];
            }                        
        }
    };
}