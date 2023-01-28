swig -lua -c++ -Iinclude cdiodeladder.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o cdiodeladder.so cdiodeladder_wrap.cxx  -lstdc++ -lm -lluajit
