    swig -lua -c++ -Iinclude WaveShaperATanSoftClip.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o WaveShaperATanSoftClip.so WaveShaperATanSoftClip_wrap.cxx -lstdc++ -lm -lluajit    
    