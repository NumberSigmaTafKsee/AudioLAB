kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude canalog_cvaonepole.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -mfma -fopenmp -pthread -shared \
-o cvaonepole.so canalog_cvaonepole_wrap.cxx  -lstdc++ -lm -lluajit
