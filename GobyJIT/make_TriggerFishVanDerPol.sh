    swig -lua -c++ -Iinclude TriggerFishVanDerPol.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o TriggerFishVanDerPol.so TriggerFishVanDerPol_wrap.cxx -lstdc++ -lm -lluajit    
    