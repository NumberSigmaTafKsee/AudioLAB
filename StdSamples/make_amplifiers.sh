swig -lua -c++ -Iinclude Amplifiers.i
gcc -fmax-errors=1 -std=c++17 -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -fopenmp -pthread -shared \
-o Amplifiers.so Amplifiers_wrap.cxx \
-lstdc++ -lm -lluajit
