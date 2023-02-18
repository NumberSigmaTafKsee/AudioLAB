#kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_amplifiers.i
gcc -fmax-errors=1 -std=c++17 -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o audiodsp_amplifiers.so $kit/audiodsp_amplifiers_wrap.cxx \
-lstdc++ -lm -lluajit
