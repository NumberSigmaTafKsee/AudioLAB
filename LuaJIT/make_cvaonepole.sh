swig -lua -c++ -Iinclude cvaonepole.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o cvaonepole.so cvaonepole_wrap.cxx  -lstdc++ -lm -lluajit
