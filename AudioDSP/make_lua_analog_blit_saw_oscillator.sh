    swig -lua -c++ -Iinclude lua_analog_blit_saw_oscillator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o lua_analog_blit_saw_oscillator.so lua_analog_blit_saw_oscillator_wrap.cxx -lstdc++ -lm -lluajit 
    