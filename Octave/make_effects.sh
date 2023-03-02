gcc -g -std=c++17 -fmax-errors=1 -O2 -march=native -mavx2 -Iinclude/DaisySP \
-Iinclude -ISynthesizer -I/usr/local/include  -I/usr/local/include/lilv-0 -I/usr/local/include -I. -I/usr/local/include/luajit-2.1 \
-o effects audio_effects.cpp AudioMidi/audiosystem.c  lib/libfv3_float.a lib/libsr2_float.a \
lib/libfv3_double.a lib/libsamplerate2.a lib/libgdither.a lib/libsr2_double.a \
lib/libstk.a lib/libGamma.a lib/libaudiofft.a lib/libfftconvolver.a lib/libDaisySP.a \
-lstdc++ -lm -lportaudio -lportmidi -lpthread -lsndfile -lluajit -lfltk -lfftw3 -lfftw3f -lpffastconv -lpfdsp -lpffft \
-llilv-0 -lsvml -lATKCore -lATKAdaptive -lATKDelay -lATKDistortion -lATKDynamic -lATKEQ -lATKIO -lATKPreamplifier \
-lATKReverberation -lATKSpecial -lATKTools -lATKUtility  -lDSPFilters
#-lkfr_dft -lkfr_io -lkfr_capi
