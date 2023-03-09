#pragma once

#include "audiodsp_fftw.hpp"

namespace AudioDSP
{
       // Code adapted from gsl/fft/factorize.c
    void factorize (const int n,
                    int *n_factors,
                    int factors[],
                    int * implemented_factors)
    {
        int nf = 0;
        int ntest = n;
        int factor;
        int i = 0;

        if (n == 0)
        {
            printf("Length n must be positive integer\n");
            return ;
        }

        if (n == 1)
        {
            factors[0] = 1;
            *n_factors = 1;
            return ;
        }

        /* deal with the implemented factors */

        while (implemented_factors[i] && ntest != 1)
        {
            factor = implemented_factors[i];
            while ((ntest % factor) == 0)
            {
                ntest = ntest / factor;
                factors[nf] = factor;
                nf++;
            }
            i++;
        }

        // Ok that's it
        if(ntest != 1)
        {
            factors[nf] = ntest;
            nf++;
        }

        /* check that the factorization is correct */
        {
            int product = 1;

            for (i = 0; i < nf; i++)
            {
                product *= factors[i];
            }

            if (product != n)
            {
                printf("factorization failed");
            }
        }

        *n_factors = nf;
    }



    bool is_optimal(int n, int * implemented_factors)
    {
        // We check that n is not a multiple of 4*4*4*2
        if(n % 4*4*4*2 == 0)
            return false;

        int nf=0;
        int factors[64];
        int i = 0;
        factorize(n, &nf, factors,implemented_factors);

        // We just have to check if the last factor belongs to GSL_FACTORS
        while(implemented_factors[i])
        {
            if(factors[nf-1] == implemented_factors[i])
                return true;
            ++i;
        }
        return false;
    }

    int find_closest_factor(int n)
    {
        int j;
        int FFTW_FACTORS[7] = {13,11,7,5,3,2,0}; // end with zero to detect the end of the array
        if(is_optimal(n,FFTW_FACTORS))
            return n;
        else
        {
            j = n+1;
            while(!is_optimal(j,FFTW_FACTORS))
                ++j;
            return j;
        }
    }

    
    
    struct FFTConvolutionDouble
    {
        FFTPlanRealDouble fftPlan;
        size_t  length,fftSize,blockSize;
 
        std::vector<std::complex<double>> t1,t2,tempC;
        std::vector<double> H,i1,ola,temp;
        
        FFTConvolutionDouble(size_t len, double * h, size_t blocks)
        {           
            length  = len;            
            blockSize = blocks;
            fftSize = find_closest_factor(len + blocks -1);
            
            fftPlan.init(fftSize);            
            
            t1.resize(fftSize);
            t2.resize(fftSize);
            tempC.resize(fftSize);

            H.resize(fftSize);
            memset(H.data(),0,H.size()*sizeof(double));
            memcpy(H.data(),h,len*sizeof(double));     

            i1.resize(fftSize);
            memset(i1.data(),0,i1.size()*sizeof(double));

            ola.resize(fftSize);
            memset(ola.data(),0,ola.size()*sizeof(double));
            
            temp.resize(fftSize);            
        }

        void ProcessBlock(size_t n, double * in, double * out)
        {
            memcpy(i1.data(),in,n*sizeof(double));

            
            fft(fftPlan,H.data(),t1.data(),false);
            fft(fftPlan,i1.data(),t2.data(),false);
                        
            tempC[0] = std::complex<double>(0,0);
            tempC[t2.size()/2-1] = std::complex<double>(0,0);

            for(size_t i = 1; i < t1.size()/2-1; i++)
            {
                tempC[i] = t1[i] * t2[i];
            }
              
            ifft(fftPlan,tempC.data(),temp.data(),true);

            for(size_t i = 0; i < n; i++)
            {
                out[i] = (temp[i] + ola[i]);                
                ola[i] = ola[i+n];                
            }            
            for(size_t i = n; i < temp.size(); i++)
                ola[i] = temp[i];            
        }
    };

    struct FFTConvolutionFloat
    {
        FFTPlanRealFloat fftPlan;
        size_t  length,fftSize,blockSize;
 
        std::vector<std::complex<float>> t1,t2,tempC;
        std::vector<float> H,i1,ola,temp;
        
        FFTConvolutionFloat(size_t len, float * h, size_t blocks)
        {           
            length  = len;            
            blockSize = blocks;
            fftSize = find_closest_factor(len + blocks -1);
            
            fftPlan.init(fftSize);            
            
            t1.resize(fftSize);
            t2.resize(fftSize);
            tempC.resize(fftSize);

            H.resize(fftSize);
            memset(H.data(),0,H.size()*sizeof(float));
            memcpy(H.data(),h,len*sizeof(float));     

            i1.resize(fftSize);
            memset(i1.data(),0,i1.size()*sizeof(float));

            ola.resize(fftSize);
            memset(ola.data(),0,ola.size()*sizeof(float));
            
            temp.resize(fftSize);            
        }

        void ProcessBlock(size_t n, float * in, float * out)
        {
            memcpy(i1.data(),in,n*sizeof(float));

            
            fft(fftPlan,H.data(),t1.data(),false);
            fft(fftPlan,i1.data(),t2.data(),false);
                        
            tempC[0] = std::complex<double>(0,0);
            tempC[t2.size()/2-1] = std::complex<float>(0,0);

            for(size_t i = 1; i < t1.size()/2-1; i++)
            {
                tempC[i] = t1[i] * t2[i];
            }
              
            ifft(fftPlan,tempC.data(),temp.data(),true);

            for(size_t i = 0; i < n; i++)
            {
                out[i] = (temp[i] + ola[i]);                
                ola[i] = ola[i+n];                
            }            
            for(size_t i = n; i < temp.size(); i++)
                ola[i] = temp[i];            
        }
    };
}