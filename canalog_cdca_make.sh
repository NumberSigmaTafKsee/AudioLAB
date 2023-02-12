kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude canalog_cdca.i
gcc -fmax-errors=1 -std=c++17 -Iinclude -I.-O2 -fPIC -march=native -mavx2 -mfma -fopenmp -pthread -shared -o canalog_cdca.so \
canalog_cdca_wrap.cxx  -lstdc++ -lm -lluajit
