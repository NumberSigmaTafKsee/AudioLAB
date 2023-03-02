    swig -lua -c++ -Iinclude AudioEffectPitchShifter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectPitchShifter.so AudioEffectPitchShifter_wrap.cxx -lstdc++ -lm -lluajit    
    