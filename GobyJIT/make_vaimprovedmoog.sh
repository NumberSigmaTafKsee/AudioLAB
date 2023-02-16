swig -lua -c++ vaimprovedmoog.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o vaimprovedmoog.so vaimprovedmoog_wrap.cxx -lstdc++ -lm -lluajit
