kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_canalog.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o canalog.so \
%kit/audiodsp_canalog_wrap.cxx -lstdc++ -lm -lluajit
