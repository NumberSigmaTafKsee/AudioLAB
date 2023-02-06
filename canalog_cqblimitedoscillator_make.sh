kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_cqblimitedoscillator.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared \
-o cqblimitedoscillator.so $kit/audiodsp_cqblimitedoscillator_wrap.cxx  -lstdc++ -lm -lluajit
