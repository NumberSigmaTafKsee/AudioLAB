swig -lua -c++ jack.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o jack.so jack_wrap.cxx jack.cpp -lstdc++ -lm -lluajit -lsndfile
