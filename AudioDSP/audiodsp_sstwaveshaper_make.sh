swig -lua -c++ -Iinclude Kits/AudioDSP/audiodsp_sstwaveshaper.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o sstwaveshaper.so \
Kits/AudioDSP/audiodsp_sstwaveshaper_wrap.cxx  -lstdc++ -lm -lluajit
