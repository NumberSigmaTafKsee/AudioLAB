swig -lua -c++ -Iinclude -D__UNIX_JACK__ -D__LINUX_PULSE__ -DRAWWAVE_PATH=Data/rawwaves Stk.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I/usr/local/include/lilv-0 \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o Stk.so Stk_wrap.cxx lib/libstk.a \
-lstdc++ -lm -lluajit
