swig -lua -c++ -Iinclude audiodsp_ck35filter.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 \
-shared -o ck35filter.so audiodsp_ck35filter_wrap.cxx  -lstdc++ -lm -lluajit
