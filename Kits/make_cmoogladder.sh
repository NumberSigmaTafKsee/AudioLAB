swig -lua -c++ -Iinclude audiodsp_cmoogladder.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared \
-o cmoogladder.so audiodsp_cmoogladder_wrap.cxx  -lstdc++ -lm -lluajit
