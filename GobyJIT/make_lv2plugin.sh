swig -lua -c++ lv2plugin.i
gcc -g -std=c++17 -fmax-errors=1 -O2 -march=native -mavx2 -fPIC -shared \
-Iinclude -I/usr/local/include  -I/usr/local/include/lilv-0 -I/usr/local/include -I. -I/usr/local/include/luajit-2.1 \
-o lv2plugin.so lv2plugin_wrap.cxx -lstdc++ -lm -lpthread -lsndfile -lluajit -llilv-0 