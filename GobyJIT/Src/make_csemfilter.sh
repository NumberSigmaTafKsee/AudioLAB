swig -lua -c++ -Iinclude csemfilter.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared \
-o csemfilter.so csemfilter_wrap.cxx  -lstdc++ -lm -lluajit
