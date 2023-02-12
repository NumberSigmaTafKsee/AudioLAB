kit="Kits/AudioDSP"
swig -lua -c++ -Iinclude canalog_cdiodeladder.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -mfma -fopenmp -pthread -shared \
-o cdiodeladder.so canalog_cdiodeladder_wrap.cxx  -lstdc++ -lm -lluajit
