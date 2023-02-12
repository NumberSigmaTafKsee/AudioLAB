swig -lua -c++ -Iinclude analog_diode.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o analog_diode.so analog_diode_wrap.cxx -lstdc++ -lm -lluajit
