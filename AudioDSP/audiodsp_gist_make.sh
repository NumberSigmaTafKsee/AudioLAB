kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude/Gist -Iinclude/Gist/src $kit/audiodsp_gist.i
gcc -Iinclude/Gist -Iinclude/Gist/src -O2 -fPIC -march=native -mavx2 -shared -o gist.so \
$kit/audiodsp_gist_wrap.cxx lib/libGist.a -lstdc++ -lm -lluajit -L. -lfftw3
