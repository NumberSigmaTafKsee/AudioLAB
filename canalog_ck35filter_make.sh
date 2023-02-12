kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude canalog_ck35filter.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -mfma -fopenmp -pthread \
-shared -o ck35filter.so canalog_ck35filter_wrap.cxx  -lstdc++ -lm -lluajit
