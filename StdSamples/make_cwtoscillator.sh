swig -lua -c++ -Iinclude audiodsp_cwtoscillator.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared \
-o cwtoscillator.so audiodsp_cwtoscillator_wrap.cxx  -lstdc++ -lm -lluajit
