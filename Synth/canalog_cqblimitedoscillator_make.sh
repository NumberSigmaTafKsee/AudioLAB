kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude canalog_cqblimitedoscillator.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -mfma -fopenmp -pthread -shared \
-o cqblimitedoscillator.so canalog_cqblimitedoscillator_wrap.cxx  -lstdc++ -lm -lluajit
