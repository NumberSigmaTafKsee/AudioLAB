    swig -lua -c++ -Iinclude FxDSPSpectrumAnalyzer.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPSpectrumAnalyzer.so FxDSPSpectrumAnalyzer_wrap.cxx -lstdc++ -lm -lluajit    
    