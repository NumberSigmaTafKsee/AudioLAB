%module kissfft
%{ 
#include <kiss_fft.h>
#include <cassert>
#include <complex>

using namespace std;
%}

//%include <kiss_fft.h>
%include "std_math.i"
%include "stdint.i"
%include "std_vector.i"

using namespace std;

%native(table_to_vector) int TableToVector(lua_State *L);
%native(vector_to_table) int VectorToTable(lua_State *L);
%{ 
int TableToVector(lua_State *L)
{
    size_t l = lua_objlen(L,-1);        
    std::vector<float> r(l);
    for(size_t i = 0; i < l; i++)
    {
        lua_rawgeti(L,-1,i+1);
        float v = lua_tonumber(L,-1);        
        r[i] = v;        
        lua_pop(L,1);
    }
    swig_type_info * info = swig_types[11];
    std::vector<float> * p = new std::vector<float>(l);
    for(size_t i = 0; i < l; i++)
        (*p)[i] = r[i];
    SWIG_NewPointerObj(L,p,info,1); 
    return 1;
}
int VectorToTable(lua_State * L)
{
    std::vector<float> *v;
    swig_type_info * info = swig_types[11];
    SWIG_ConvertPtr(L,-1,(void**)&v,info,1);
    lua_createtable(L,v->size(),0);
    for(size_t i = 0; i < v->size(); i++)
    {
        lua_pushnumber(L,(*v)[i]);
        lua_rawseti(L,-2,i+1);    
    }
    return 1;
}
%}

%inline 
%{
    struct KissFFT
    {
        kiss_fft_cfg  config;
        kiss_fft_cfg  iconfig;
        
        KissFFT(size_t nfft)
        {
            config = kiss_fft_alloc(nfft,0,NULL,NULL);
            assert(config != NULL);
            iconfig = kiss_fft_alloc(nfft,1,NULL,NULL);
            assert(iconfig != NULL);        
        }
        ~KissFFT() 
        {
            if(config) free(config);
            if(iconfig) free(iconfig);
        }
        
        void stride(vector<kiss_fft_cpx> & in, vector<kiss_fft_cpx> & out, size_t stride)
        {
            out.resize(in.size());
            kiss_fft_stride(config,in.data(), out.data(), stride);
        }

        void forward(std::complex<float>* in, 
                                std::complex<float>* out)
        {                 
            kiss_fft(config, (const kiss_fft_cpx*)in, (kiss_fft_cpx*)out);     
        }
        void inverse(std::complex<float>* in, 
                                std::complex<float>* out)
        {                    
            kiss_fft(iconfig, (const kiss_fft_cpx*)in, (kiss_fft_cpx*)out);        
        }
    };
%}