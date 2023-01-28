swig -lua -c++ -Iinclude -Iinclude Kits/AudioDSP/audiodsp_stk.i 
gcc -D__UNIX_JACK__ -D__LINUX_PULSE__ -DRAWWAVE_PATH=Data/rawwaves \
-fmax-errors=1 -std=c++17 -I. -Iinclude -IKits/AudioDSP \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o Stk.so Kits/AudioDSP/audiodsp_stk_wrap.cxx lib/libstk.a \
-lstdc++ -lm -lluajit
