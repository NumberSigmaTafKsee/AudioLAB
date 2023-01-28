kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_fxdistortion.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o fxdistortion.so \
$kit/audiodsp_fxdistortion_wrap.cxx -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
