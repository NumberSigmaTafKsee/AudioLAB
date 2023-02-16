#pragma once

namespace KfrDSP1
{
///////////////////////////////////////////////////////////////
// Convolution
///////////////////////////////////////////////////////////////

    template<typename T> using ConvolveFilter = DSP::ConvolveFilter<T>;
    
    template<typename T>
    struct ConvolutionFilter
    {
        kfr::convolve_filter<T> *filter;        
        kfr::univector<T>        temp,out,ola;
        size_t h_size;
        size_t block_size;        
        ConvolutionFilter(kfr::univector<T> &impulse,size_t block_size_=1024)
        {
            filter = new kfr::convolve_filter<T>(impulse,block_size);                   
            temp.resize(block_size);
            out.resize(block_size);
            ola.resize(block_size);
            memset(ola.data(),0,block_size*sizeof(T));
            h_size = impulse.size();
            block_size = block_size_;
        }
        ~ConvolutionFilter() {
            if(filter) delete filter;
        }
        void ProcessBlock(size_t n, T * input, T * output) {
            memcpy(temp.data(),input,n*sizeof(T));
            filter->apply(temp,out);
            for(size_t i = 0; i < n; i++) {
                output[i] = out[i];
            }
            if(n < block_size) {
                for(size_t i = 0; i < n; i++) output[i] += ola[i];
                for(size_t i = n; i < block_size; i++)                                    
                    ola[i-n] = ola[i];
            }
        }
        kfr::univector<T> Process(kfr::univector<T> & m) 
        {
            kfr::univector<T> r(m.size()+h_size-1);                        
            filter->apply(r,m);            
            return r;
        }
    };

    template<typename T>
    struct StereoConvolutionFilter
    {
        ConvolutionFilter<T> *filter[2];        
        kfr::univector<T>        temp[2];
        kfr::univector<T>        out[2];
        StereoConvolutionFilter(kfr::univector<T> &impulseL, kfr::univector<T> &impulseR,size_t block_size=1024)
        {
            filter[0] = new ConvolutionFilter<T>(impulseL,block_size);  
            filter[1] = new ConvolutionFilter<T>(impulseR,block_size);                   
            temp[0].resize(block_size);
            temp[1].resize(block_size);
            out[0].resize(block_size);
            out[1].resize(block_size);
        }
        ~StereoConvolutionFilter() {
            if(filter[0]) delete filter[0];
            if(filter[1]) delete filter[1];
        }
        void ProcessBlock(size_t n, T ** in, T ** out) {
            filter[0]->ProcessBlock(n,in[0],out[0]);
            filter[1]->ProcessBlock(n,in[1],out[1]);
        }        
        kfr::univector2d<T> Process(kfr::univector2d<T> & m) 
        {
            kfr::univector2d<T> r(2);            
            
            r[0] = filter[0]->Process(m[0]);
            r[1] = filter[1]->Process(m[1]);
            return r;
        }
    };


    //////////////////////////////////////////////////////////////
    // Convolve/Correlation
    ///////////////////////////////////////////////////////////////

    template<typename T>
    kfr::univector<T> convolve(kfr::univector<T> a, kfr::univector<T> b)
    {                
        return kfr::convolve(a,b);        
    }
    template<typename T>
    kfr::univector<T> correlate(kfr::univector<T> a, kfr::univector<T> b)
    {                
        return kfr::correlate(a,b);        
    }
    template<typename T>
    kfr::univector<T> autocorrelate(kfr::univector<T> a)
    {                
        return kfr::autocorrelate(a);        
    }
}        