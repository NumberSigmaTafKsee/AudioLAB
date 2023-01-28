swig -octave -c++ -Iinclude audio_dc_block_filter.i
mkoctfile -Wfatal-errors -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -mavx2 -mfma -march=native \
-o audio_dc_block_filter audio_dc_block_filter_wrap.cxx -lstdc++ -lm -lluajit
