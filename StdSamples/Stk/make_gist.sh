swig -lua -c++ -IGist/Gist/src gist.i
gcc -I. -IGist/Gist/src -O2 -fPIC -march=native -mavx2 -shared -o gist.so gist_wrap.cxx lib/libGist.a -lstdc++ -lm -lluajit -L. -lfftw3
