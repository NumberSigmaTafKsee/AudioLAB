swig -lua -c++ -Iinclude audiodsp_analog.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o audiodsp_analog.so audiodsp_analog_wrap.cxx \
-lstdc++ -lm -lluajit -lutil -ldl -lm -lpthread
