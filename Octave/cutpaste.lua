template = [[
~module %s
~{
typedef float DspFloatType;
#include "SoundObject.hpp"
#include "%s"
~}
typedef float DspFloatType;
~include "stdint.i"
~include "std_math.i"
~include "std_vector.i"
~include "std_list.i"
~include "std_map.i"

typedef float DspFloatType;
~include "SoundObject.hpp"
~include "%s"

~template(float_vector) std::vector<float>;
~template(double_vector) std::vector<double>;

~template(complex_float_vector) std::vector<std::complex<float>>;
~template(complex_double_vector) std::vector<std::complex<double>>;

]]

make = [[
    swig -lua -c++ -I../include -I../include %s.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o %s.so %s_wrap.cxx -lstdc++ -lm -lluajit 
    ]]
    
function make_files(dir,hdr)
    x = string.format(template,hdr, "FX/" .. hdr .. ".hpp", "FX/" .. hdr .. ".hpp")
    x = x:gsub("~","%%")    
    f = io.open(hdr .. ".swg",'w')
    f:write(x)
    f:close()

    x = string.format(make,hdr,hdr,hdr)
    f = io.open("make_" .. hdr .. ".sh",'w')
    f:write(x)
    f:close()
    --os.execute("sh make_" .. hdr .. ".sh")
end

require("lfs")

for file in lfs.dir(arg[1])  do	
	print(file)    
	if(file:find("hpp") ~= 0) then		
		make_files(arg[1],file:gsub(".hpp",""))    
	end
end
