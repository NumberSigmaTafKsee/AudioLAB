swig -Isrc -lua -c++ -Iinclude -Ikfrcore/include/Kfr src/SoundWave.i
gcc -std=c++17 -Isrc -Iinclude -Iinclude/samples -Llib -DAUDIOFFT_FFTW3 -fmax-errors=1 -I/usr/local/include/kissfft -O2 -march=native -mavx2 -msse -fPIC -shared -o soundwave.so src/SoundWave_wrap.cxx lib/libSoundWave.a -lstdc++ -lm -lluajit -L../lib -lsndfile -lkissfft-float -lsamplerate
