swig -lua -c++ -Iinclude Envelopes.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o Envelopes.so Envelopes_wrap.cxx  \
-lstdc++ -lm -lluajit
