swig -lua -c++ -Iinclude/Gist src/gist.i
gcc -I. -Iinclude/Gist -O2 -fPIC -march=native -mavx2 -shared -o gist.so src/gist_wrap.cxx lib/libGist.a -lstdc++ -lm -lluajit -L. -lfftw3
