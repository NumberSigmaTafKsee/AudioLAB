swig -lua -c++ -Iinclude audiodsp_cqblimitedoscillator.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared \
-o cqblimitedoscillator.so audiodsp_cqblimitedoscillator_wrap.cxx  -lstdc++ -lm -lluajit
