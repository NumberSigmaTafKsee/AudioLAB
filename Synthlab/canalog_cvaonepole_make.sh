kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_cvaonepole.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared \
-o cvaonepole.so $kit/audiodsp_cvaonepole_wrap.cxx  -lstdc++ -lm -lluajit
