    swig -lua -c++ -Iinclude audiodsp_fft_wavetable_generator.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_fft_wavetable_generator.so audiodsp_fft_wavetable_generator_wrap.cxx -lstdc++ -lm -lluajit 
    