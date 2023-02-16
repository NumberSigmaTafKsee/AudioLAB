swig -lua -c++ -Iinclude Analog.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o Analog.so Analog_wrap.cxx \
-lstdc++ -lm -lluajit -lutil -ldl -lm -lpthread
