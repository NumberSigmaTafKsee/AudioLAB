swig -lua -c++ mkl_array.i
gcc -fmax-errors=1 -std=c++17  -DMKL_ILP64  -m64  -I"${MKLROOT}/include" -O2 -march=native -mavx2 -fPIC -shared -o mkl_array.so  mkl_array_wrap.cxx -lstdc++ -lm -lluajit  -L/opt/intel/oneapi/mkl/2022.0.2/lib/intel64 -Wl,--no-as-needed  -lmkl_rt -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread -lm -ldl -llapack
