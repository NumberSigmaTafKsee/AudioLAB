kit="Kits/AudioDSP"
kit="."
swig -lua -c++ $kit/audiodsp_fftconvolver.i
gcc -DAUDIOFFT_FFTW3 -O2 -fPIC -march=native -mavx2 -shared -o fftconvolver.so \
$kit/audiodsp_fftconvolver_wrap.cxx lib/libfftconvolver.a -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
