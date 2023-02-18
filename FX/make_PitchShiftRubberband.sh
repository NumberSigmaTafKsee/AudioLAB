    swig -lua -c++ -I../include -I../include PitchShiftRubberband.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o PitchShiftRubberband.so PitchShiftRubberband_wrap.cxx -lstdc++ -lm -lluajit 
    