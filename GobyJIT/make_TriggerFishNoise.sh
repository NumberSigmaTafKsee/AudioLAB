    swig -lua -c++ -Iinclude TriggerFishNoise.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o TriggerFishNoise.so TriggerFishNoise_wrap.cxx -lstdc++ -lm -lluajit    
    