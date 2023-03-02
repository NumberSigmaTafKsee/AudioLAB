clear
gcc -fmax-errors=1 -Iinclude -Iinclude/SimpleEigen -O2 -o xor cranium.cpp -lstdc++ -lm -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_thread -lmkl_rt -lmkl_core -lmkl_intel_lp64  -liomp5 -lpthread -lm -ldl
