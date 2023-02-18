    swig -lua -c++ -Iinclude audiodsp_rubberband_pitch_shifter.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_rubberband_pitch_shifter.so audiodsp_rubberband_pitch_shifter_wrap.cxx -lstdc++ -lm -lluajit 
    