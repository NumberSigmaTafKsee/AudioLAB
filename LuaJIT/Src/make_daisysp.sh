swig -lua -c++ -Iinclude -Iinclude/DaisySP DaisySP.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude/DaisySP -Iinclude -I/usr/local/include/lilv-0 \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o DaisySP.so DaisySP_wrap.cxx lib/libDaisySP.a \
-lstdc++ -lm -lluajit 
