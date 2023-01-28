swig -lua -Iinclude Kits/AudioDSP/audiodsp_xtract.i
gcc -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o xtract.so Kits/AudioDSP/audiodsp_xtract_wrap.c -lm -lluajit
