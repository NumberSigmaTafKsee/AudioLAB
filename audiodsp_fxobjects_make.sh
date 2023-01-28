kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_fxobjects.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o fxobjects.so \
$kit/audiodsp_fxobjects_wrap.cxx -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
