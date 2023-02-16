swig -lua -c++ vakorg35hpf.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o vakorg35hpf.so vakorg35hpf_wrap.cxx -lstdc++ -lm -lluajit
