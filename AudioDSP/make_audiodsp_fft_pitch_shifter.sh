    swig -lua -c++ -Iinclude audiodsp_fft_pitch_shifter.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_fft_pitch_shifter.so audiodsp_fft_pitch_shifter_wrap.cxx -lstdc++ -lm -lluajit 
    