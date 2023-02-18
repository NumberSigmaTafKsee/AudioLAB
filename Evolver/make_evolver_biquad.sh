gcc -std=c++17 -O3 -Iinclude -march=native -mavx2 -mfma -o evolver_biquad evolver_biquad.cpp ga.c mt19937.c -lstdc++ -lm -lfftw3 -lfftw3f
