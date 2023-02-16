swig -lua -c++ -IGist/Gist/src audiodsp_gist.i
gcc -I. -IGist/Gist/src -O2 -fPIC -march=native -mavx2 -shared -o gist.so audiodsp_gist_wrap.cxx lib/libGist.a -lstdc++ -lm -lluajit -L. -lfftw3
