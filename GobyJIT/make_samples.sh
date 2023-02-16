swig -lua -c++ -Iinclude/samples -Iinclude/Kfr src/samples.i
gcc -std=c++17 -fmax-errors=1 -I. -Iinclude -Iinclude/Kfr -Iinclude/samples -O2 -fPIC -shared -o samples.so src/samples_wrap.cxx -lstdc++ -lm -lluajit
