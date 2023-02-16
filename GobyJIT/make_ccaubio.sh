swig -lua -c++ -Iaubio/src ccaubio.i
gcc -Iaubio/src -O2 -fPIC -march=native -mavx2 -shared -o ccaubio.so ccaubio_wrap.cxx ../lib/libaubio.a -lstdc++ -lm -lluajit -lsndfile -lsamplerate -lavformat -lrubberband
