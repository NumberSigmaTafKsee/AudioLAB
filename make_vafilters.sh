swig -lua -c++ -Iinclude Kits/AudioDSP/audiodsp_vafilters.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o vafilters.so \
Kits/AudioDSP/audiodsp_vafilters_wrap.cxx -lstdc++ -lm -lluajit
