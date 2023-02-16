swig -lua -c++ -Iaubio aubio.i
gcc -Iaubio -fpermissive -Iaubio/src -O2 -fPIC -march=native -mavx2 -shared -o aubio.so aubio_wrap.cxx lib/libaubio.a -lstdc++ -lm -lluajit -lsndfile -lsamplerate -lavformat -lrubberband
