swig -lua -c++ -I../include -I.. amplifiersfuncs.i
gcc -fmax-errors=1 -I.. -I../include -O2 -fPIC -std=c++17 -march=native -mavx2 -shared -o amplifiersfuncs.so amplifiersfuncs_wrap.cxx -lstdc++ -lm -lluajit
