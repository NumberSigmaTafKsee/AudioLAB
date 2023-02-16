swig -lua -c++ -Iccaubio/aubio ccaubio.i
gcc -Iccaubio/aubio -O2 -fPIC -march=native -mavx2 -shared -o ccaubio.so ccaubio_wrap.cxx lib/libaubio.a -lstdc++ -lm -lluajit -lsndfile -lsamplerate -lavformat -lrubberband
