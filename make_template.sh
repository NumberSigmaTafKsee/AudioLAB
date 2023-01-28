gcc -g -std=c++17 -fmax-errors=1 -O2 -march=native -mavx2 -Iinclude/DaisySP -I/usr/local/include/octave-7.3.0 \
-Iinclude -ISynthesizer -I/usr/local/include  -I/usr/local/include -I. -I/usr/local/include/luajit-2.1 \
-o template audio_template.cpp AudioMidi/audiosystem.c  lib/libfv3_float.a lib/libsr2_float.a \
lib/libfv3_double.a lib/libsamplerate2.a lib/libgdither.a lib/libsr2_double.a \
lib/libstk.a lib/libGamma.a lib/libaudiofft.a lib/libfftconvolver.a lib/libDaisySP.a \
-lstdc++ -lm -lportaudio -lportmidi -lpthread -lsndfile -lluajit -lfltk -lfftw3 -lfftw3f -lpffastconv -lpfdsp -lpffft \
-llilv-0 -lsvml -lATKCore -lATKAdaptive -lATKDelay -lATKDistortion -lATKDynamic -lATKEQ -lATKIO -lATKPreamplifier \
-lATKReverberation -lATKSpecial -lATKTools -lATKUtility  -lDSPFilters -L/usr/local/lib/octave/7.3.0 -loctave -loctinterp
#-lkfr_dft -lkfr_io -lkfr_capi
