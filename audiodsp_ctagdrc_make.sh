swig -lua -c++ -Iinclude audiodsp_ctagdrc.i
gcc -fmax-errors=1 -std=c++17 -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o audiodsp_ctagdrc.so audiodsp_ctagdrc_wrap.cxx \
-lstdc++ -lm -lluajit
