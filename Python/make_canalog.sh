swig -lua -c++ -Iinclude -ICAnalog audiodsp_canalog.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -ICAnalog -I.-O2 -fPIC -march=native -mavx2 -shared -o canalog.so audiodsp_canalog_wrap.cxx -lstdc++ -lm -lluajit
