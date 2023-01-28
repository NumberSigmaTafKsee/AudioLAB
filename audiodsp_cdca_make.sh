kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_cdca.i
gcc -fmax-errors=1 -std=c++17 -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o cdca.so \
Kits/AudioDSP/audiodsp_cdca_wrap.cxx  -lstdc++ -lm -lluajit
