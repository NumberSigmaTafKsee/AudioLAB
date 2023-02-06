kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_audiotk.i
gcc -fmax-errors=1 -std=c++17 -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o AudioTK.so $kit/audiodsp_audiotk_wrap.cxx \
-lstdc++ -lm -lluajit -lATKCore_static -lATKSpecial_static -lATKAdaptive_static -lATKDelay_static -lATKDistortion_static -lATKDynamic_static \
-lATKEQ_static -lATKIO_static -lATKPreamplifier_static -lATKReverberation_static -lATKTools_static -lATKUtility_static -lfftw3 -lfftw3f -lsndfile -lsamplerate
