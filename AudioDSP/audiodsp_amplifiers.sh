swig -lua -c++ -Iinclude audiodsp_amplifiers.i
gcc -fmax-errors=1 -std=c++17 -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -fopenmp -pthread -shared \
-o audiodsp_amplifiers.so audiodsp_amplifiers_wrap.cxx \
-lstdc++ -lm -lluajit
