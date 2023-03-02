swig -lua -c++ -Iinclude src/TKorg35LPF.i
gcc -fmax-errors=1 -Iinclude -O2 -march=native -mavx2 -fPIC -shared -o TKorg35LPF.so src/TKorg35LPF_wrap.cxx -lstdc++ -lm -lluajit
