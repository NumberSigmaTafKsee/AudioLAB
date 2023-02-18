swig -octave -c++ -Iinclude audiodsp_analog_svf.i
mkoctfile -Wfatal-errors -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -mavx2 -mfma -march=native \
-o audiodsp_analog_svf audiodsp_analog_svf_wrap.cxx -lstdc++ -lm -lluajit
