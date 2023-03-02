swig -python -c++ -Iinclude vector.i
gcc -fmax-errors=1 -Iinclude -I/usr/local/include/python3.9 -O2 -fPIC -march=native -mavx2 -shared -o _vector.so vector_wrap.cxx -lstdc++ -lm -lpython3.9
