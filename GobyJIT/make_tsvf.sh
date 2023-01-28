swig -lua -c++ -Iinclude src/TStateVariableFilter.i
gcc -Iinclude -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o TStateVariableFilter.so src/TStateVariableFilter_wrap.cxx -lstdc++ -lm -lluajit
