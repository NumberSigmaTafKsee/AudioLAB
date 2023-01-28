#gcc  -fmax-errors=1 -fopenmp -pthread -std=c++17 -I. -Iinclude -O2 -mavx2 -mfma  -march=native \
#     -o cranium cranium.cpp -lstdc++ -lm -lblas -llapacke


gcc -DUSE_MKL -fmax-errors=1 -DMKL_ILP64-m64 -mavx2 -mfma  -I"${MKLROOT}/include"-Iinclude/SimpleEigen -O2 -march=native -mavx2 -o cranium cranium.cpp -lstdc++ -L/usr/local/cuda/lib64 -lnvblas -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl -lfftw3 -lfftw3f -lsndfile  -lblas -llapacke

#gcc  -fmax-errors=1 -fopenmp -pthread -std=c++17 -I. -Iinclude -O2 -mavx2 -mfma  -march=native \
 #    -o cranium cranium.cpp -lstdc++ -lm -L/usr/local/cuda/lib64 -lnvblas -lopenblas -llapacke
     

#gcc -DUSE_MKL -fmax-errors=1 -DMKL_ILP64-m64 -mavx2 -mfma  -I"${MKLROOT}/include"-Iinclude/SimpleEigen \
#     -O2 -march=native -mavx2 -o cranium cranium.cpp -lstdc++ -L. -lblas2cuda -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl -lfftw3 -lfftw3f -lsndfile  -lblas -llapacke

