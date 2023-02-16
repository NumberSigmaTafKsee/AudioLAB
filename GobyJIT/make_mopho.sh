gcc -std=c++17 -fmax-errors=1 -O2 -march=native -mavx2 \
-IAudioMidi -Iinclude -I/usr/local/include/lilv-0 -I/usr/local/include -I. -I/usr/local/include/luajit-2.1 \
-o mopho audio_mopho.cpp AudioMidi/audiosystem.c \
lib/libfv3_float.a lib/libsamplerate2.a lib/libgdither.a lib/libsr2_float.a lib/libstk.a lib/libGamma.a lib/libaudiofft.a lib/libfftconvolver.a \
-lstdc++ -lm -lportaudio -lportmidi -lkfr_dft -lkfr_io -lpthread -lsndfile -lluajit \
-lfaustwithllvm -lrt -ldl -lLLVM-10 -lz -lcurses -lfftw3 -lfftw3f -llilv-0

