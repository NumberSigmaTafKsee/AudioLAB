gcc -std=c++17 -I../include -O3 -march=native -pthread -mavx2 -mfma -o svf evolver_svf.cpp ga.c mt19937.c -lstdc++ -lm -lfftw3 -lfftw3f -lsndfile
