swig -lua -c++ -Icpp-spectrogram/include src/cpp-spectrogram.i
gcc -Icpp-spectrogram/include -O2 -fPIC -march=native -mavx2 -shared -o cppspectrogram.so src/cpp-spectrogram_wrap.cxx lib/libcpp-spectrogram.a -lstdc++ -lm -lluajit
