%module vector
%{
#include <vector>
#include "ctype.h"
%}
%include "stdint.i"
%include "lua_fnptr.i"
%include "std_vector.i"

%include "ctype.h"

%template(lua_vector) std::vector<SWIGLUA_REF>;
%template(ctype_vector) std::vector<CType>;

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

/*
%inline %{
    
    void for_each(std::vector<SWIGLUA_REF> & v, SWIGLUA_REF fn) {
        
        for(size_t i = 0; i < v.size(); i++)
        {
            swiglua_ref_get(&fn);    
            swiglua_ref_get(&v[i]);
            if(lua_pcall(fn.L,1,0,0) != 0)
                std::cout << "for-each error: " << lua_tostring(fn.L,-1) << std::endl;            
            lua_pop(fn.L,1);
        }
    }
%}
*/
