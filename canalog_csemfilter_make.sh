kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude canalog_csemfilter.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -mfma -fopenmp -pthread -shared \
-o csemfilter.so canalog_csemfilter_wrap.cxx  -lstdc++ -lm -lluajit
