    swig -lua -c++ -Iinclude AudioEffectRoomReverb.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectRoomReverb.so AudioEffectRoomReverb_wrap.cxx -lstdc++ -lm -lluajit    
    