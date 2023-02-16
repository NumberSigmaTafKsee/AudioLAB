require("lfs")

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
~include "lua_fnptr.i"

~include "SoundObject.hpp"
~include "%s"

~template(float_vector) std::vector<float>;
~template(double_vector) std::vector<double>;

~template(complex_float_vector) std::vector<std::complex<float>>;
~template(complex_double_vector) std::vector<std::complex<double>>;

~inline ~{
    const int BufferSize = 256;
    Std::RandomMersenne noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;
~}
]]

make = [[
    swig -lua -c++ -Iinclude %s.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o %s.so %s_wrap.cxx -lstdc++ -lm -lluajit    
    ]]
    
function make_files(dir,hdr)
    x = string.format(template,hdr,dir .."/" ..hdr .. ".hpp",dir .."/" ..hdr .. ".hpp")
    x = x:gsub("~","%%")
    f = io.open(hdr .. ".i",'w')
    f:write(x)
    f:close()

    x = string.format(make,hdr,hdr,hdr)
    f = io.open("make_" .. hdr .. ".sh",'w')
    f:write(x)
    f:close()
    os.execute("sh make_" .. hdr .. ".sh")
end

for file in lfs.dir(arg[1]) do    
    --local t = string.gsub(file,".hpp","")
    --print(file,t)
    make_files(arg[1],file:gsub(".hpp",""))    
end
