swig -lua -c++ -Iinclude Kits/AudioDSP/audiodsp_sndfile.i
gcc -std=c++17 -Iinclude  -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o sndfile.so \
Kits/AudioDSP/audiodsp_sndfile_wrap.cxx -lstdc++ -lm -lluajit -lsndfile
