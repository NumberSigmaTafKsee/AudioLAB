gcc -std=c++17 -I../include -I../StdSamples -O3 -march=native -pthread -mavx2 -mfma -o evolver evolver.cpp ga.c mt19937.c -lstdc++ -lm -lfftw3 -lfftw3f -lsndfile
