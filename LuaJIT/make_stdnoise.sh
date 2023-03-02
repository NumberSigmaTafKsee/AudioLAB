swig -lua -c++ stdnoise.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o stdnoise.so stdnoise_wrap.cxx -lstdc++ -lm -lluajit
