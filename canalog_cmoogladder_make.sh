kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude canalog_cmoogladder.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -mfma -fopenmp -pthread -shared \
-o cmoogladder.so canalog_cmoogladder_wrap.cxx  -lstdc++ -lm -lluajit
