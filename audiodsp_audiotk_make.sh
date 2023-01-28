kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_audiotk.i
gcc -fmax-errors=1 -std=c++17 -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o AudioTK.so $kit/audiodsp_audiotk_wrap.cxx \
-lstdc++ -lm -lluajit
