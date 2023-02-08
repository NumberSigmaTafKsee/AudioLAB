swig -lua -c++ -Iinclude va_analog_svf.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o va_analog_svf.so \
va_analog_svf_wrap.cxx -lstdc++ -lm -lluajit
