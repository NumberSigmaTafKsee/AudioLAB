swig -octave -c++ -I/usr/local/include -Iinclude audiosystem.i
mkoctfile -fmax-errors=1 -fpermissive -I/usr/local/include -Iinclude  -O2 -march=native -mavx2 -fPIC -shared -o audiosystem audiosystem_wrap.cxx -pthread -lrt -lm -L/usr/local/lib -lportaudio -lportmidi
