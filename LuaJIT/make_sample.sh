swig -lua -c++ -Iinclude/samples src/samples.i
gcc -std=c++17 -Iinclude -Iinclude/samples -O2 -fPIC -shared -o sample.so src/samples_wrap.cxx -lstdc++ -lm -lluajit
