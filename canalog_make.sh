kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude canalog.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o canalog.so \
canalog_wrap.cxx -lstdc++ -lm -lluajit
