swig -lua -c++ -Iinclude audiodsp_envelopes.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o Envelopes.so audiodsp_envelopes_wrap.cxx  \
-lstdc++ -lm -lluajit
