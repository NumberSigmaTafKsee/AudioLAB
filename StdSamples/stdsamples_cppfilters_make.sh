kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_cppfilters.i
gcc -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o cppfilters.so \
$kit/audiodsp_cppfilters_wrap.cxx -lstdc++ -lm -lluajit
