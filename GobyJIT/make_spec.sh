swig -lua -c++ -Iinclude cpp-spectrogram.i
gcc -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o cppspectrogram.so cpp-spectrogram_wrap.cxx libcpp-spectrogram.a -lstdc++ -lm -lluajit -lsndfile -lfreeimage
