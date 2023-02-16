%module aafilter
%{
#include "aafilter.hpp"
#include <vector>
#include <map>
#include <list>
#include <string>
using namespace Filters::AAFilter;
%}

%include "std_math.i"
%include "std_common.i"
%include "std_pair.i"
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_list.i"
%include "lua_fnptr.i"
%include "carrays.i"

%inline %{
    bool operator == (const SWIGLUA_REF a, const SWIGLUA_REF b) {
        return a.L == b.L && a.ref == b.ref;
    }
    bool operator < (const SWIGLUA_REF a, const SWIGLUA_REF b) { 
        return a.ref < b.ref;
    }
%}
%template (float_vector) std::vector<float>;
%template (double_vector) std::vector<double>;
//%template (lua_vector) std::vector<SWIGLUA_REF>;
//%template (lua_map) std::map<std::string,SWIGLUA_REF>;
//%template (lua_list) std::list<SWIGLUA_REF>;

%include "aafilter.hpp"

%template(aafilter_float) Filters::AAFilter::AAFilter<float>;
%template(aafilter_double) Filters::AAFilter::AAFilter<double>;
