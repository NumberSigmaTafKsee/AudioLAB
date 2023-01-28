swig -lua -c++ -Iinclude eigen_samples.i
gcc -fmax-errors=1 -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o sv.so eigen_samples_wrap.cxx -lstdc++ -lm -lluajit
