swig -Ilibgha -lua -c++ src/libgha.i
gcc -Ilibgha/include -O2 -fPIC -march=native -mavx2 -shared -o libgha.so src/libgha_wrap.cxx lib/libgha.a -lstdc++ -lm -lluajit
