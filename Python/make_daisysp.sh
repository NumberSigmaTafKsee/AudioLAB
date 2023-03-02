swig -lua -c++ -Iinclude -Iinclude/DaisySP audiodsp_daisysp.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude/DaisySP -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o DaisySP.so audiodsp_daisysp_wrap.cxx lib/libDaisySP.a \
-lstdc++ -lm -lluajit 
