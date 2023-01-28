#pragma once

namespace Casino::IPP
{
    template<typename T>
    void RandUniformInit(void * state, T low, T high, int seed)
    {
        assert(1==0);
    }
    template<>
    void RandUniformInit<Ipp32f>(void * state, Ipp32f low, Ipp32f high, int seed)
    {
        IppStatus status = ippsRandUniformInit_32f((IppsRandUniState_32f*)state,low,high,seed);
        checkStatus(status);
    }
    template<typename T>
    void RandUniform(T * array, int len, void * state)
    {
        assert(1==0);
    }
    template<>
    void RandUniform<Ipp32f>(Ipp32f * array, int len, void * state)
    {
        IppStatus status = ippsRandUniform_32f(array,len,(IppsRandUniState_32f*)state);
        checkStatus(status);
    }
    template<typename T>
    void RandUniformGetSize(int * p)
    {
        IppDataType dType = GetDataType<T>();
        IppStatus status;
        if(dType ==ipp32f) status = ippsRandUniformGetSize_32f(p);
        else if(dType == ipp64f) status = ippsRandUniformGetSize_64f(p);        
        checkStatus(status);
    }
    template<typename T>
    struct RandomUniform
    {
        Ipp8u * state;
        RandomUniform(T high, T low, int seed=-1)
        {
            if(seed == -1) seed = time(NULL);            
            int size=0;
            RandUniformGetSize<T>(&size);
            if(size==0) throw std::runtime_error("RandomUniform size is 0");
            state = Malloc<Ipp8u>(size);
            RandUniformInit<T>(state,high,low,seed);    
        }
        ~RandomUniform() {
            if(state) Free(state);
        }
        void fill(T* array, int len) { 
            RandUniform<T>(array,len,state);
        }
        
    };
}