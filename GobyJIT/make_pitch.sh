swig -lua -c++ pitch_detection.i
gcc -fopenmp -O2 -march=native -mavx2 -fPIC -shared -opitch_detection.so pitch_detection_wrap.cxx -lstdc++ -lluajit -lm -L. -lpitch_detection -pthread -lmlpack -lffts -larmadillo
