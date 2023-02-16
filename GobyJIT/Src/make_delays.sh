swig -lua -c++ -Iinclude Delay.i
gcc -Iinclude -fmax-errors=1 -std=c++17 -I. \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o Delay.so Delay_wrap.cxx  \
-lstdc++ -lm -lluajit
