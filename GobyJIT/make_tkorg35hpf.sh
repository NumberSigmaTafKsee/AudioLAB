swig -lua -c++ -Iinclude src/TKorg35HPF.i
gcc -Iinclude -fmax-errors=1 -Iinclude/Filters -O2 -march=native -mavx2 -fPIC -shared -o TKorg35HPF.so src/TKorg35HPF_wrap.cxx -lstdc++ -lm -lluajit
