swig -lua -c++ src/stdlite.i
gcc -fpermissive -O2 -march=native -mavx2 -fPIC -shared -o stdlite.so src/stdlite_wrap.cxx -lstdc++ -lm -lluajit
