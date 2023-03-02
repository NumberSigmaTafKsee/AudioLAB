swig -lua -c++ -Iinclude src/sndfile.i
gcc -std=c++17 -Iinclude  -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o sndfile.so src/sndfile_wrap.cxx -lstdc++ -lm -lluajit -lsndfile
