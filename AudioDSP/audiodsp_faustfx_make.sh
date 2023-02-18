kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_faustfx.i
gcc -g -std=c++17 -fmax-errors=1 -O2 -march=native -mavx2 -fPIC -shared \
-Iinclude -I/usr/local/include  -I/usr/local/include/lilv-0 -I/usr/local/include -I. -I/usr/local/include/luajit-2.1 \
-o faustfx.so $kit/audiodsp_faustfx_wrap.cxx -lstdc++ -lm -lsndfile -lluajit -lfaustwithllvm -lpthread -lrt -ldl -lLLVM-10 -lz -lcurses 
