swig -lua -c++ -Iinclude vadiodeladderfilter2.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o vadiodeladderfilter2.so vadiodeladderfilter2_wrap.cxx -lstdc++ -lm -lluajit
