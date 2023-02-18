kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude fxdistortion.i
gcc -fmax-errors=1 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o fxdistortion.so \
fxdistortion_wrap.cxx -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
