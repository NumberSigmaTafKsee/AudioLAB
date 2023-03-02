swig -DPYTHON -python -c++ -Iinclude analog_blit_oscillators.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I/usr/local/include/python3.9 -I. \
-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o _PyAnalogBlitOscillators.so \
analog_blit_oscillators_wrap.cxx -lstdc++ -lm -L/usr/local/lib -lpython3.9 -lutil -lrt -ldl
mv analog_blit_oscillators.py PyAnalogBlitOscillators.py