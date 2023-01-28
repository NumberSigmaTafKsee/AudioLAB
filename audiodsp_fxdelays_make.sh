kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_fxdelays.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o fxdelays.so \
$kit/audiodsp_fxdelays_wrap.cxx -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
