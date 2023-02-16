%module matplotlib
%{
#define WITHOUT_NUMPY    
#include "matplotlib.hpp"

namespace plt = matplotlibcpp;
%}

%include "std_math.i"
%include "std_string.i"
%include "std_vector.i"


%inline %{

    void plot(size_t n, float * y,const std::string& name="")
    {
        std::vector<float> v(n);
        plt::detail::_interpreter::get();
        memcpy(v.data(),y,n*sizeof(float));
        plt::named_plot("sin",v);            
        plt::pause(5);
    }
%}