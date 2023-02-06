swig -lua -c++ cranium.i
gcc -O2 -fPIC -march=native -mavx2 -mfma -shared -o cranium.so cranium_wrap.cxx -lstdc++ -lm -lluajit
