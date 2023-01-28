swig -lua -c++ distortionfunctions.i
gcc -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o distortionfunctions.so distortionfunctions_wrap.cxx -lstdc++ -lm -lluajit
