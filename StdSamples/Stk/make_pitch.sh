swig -lua -c++ -Iinclude pitch_detection/pitch_detection.i
gcc -Iinclude -fopenmp -O2 -march=native -mavx2 -fPIC -shared -opitch_detection.so pitch_detection/pitch_detection_wrap.cxx -Lbin -lpitch_detection -lstdc++ -lluajit -lm -pthread -lmlpack -lffts -larmadillo
