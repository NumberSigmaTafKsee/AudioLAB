swig -lua -c++ -Iinclude audio_function_generator.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o audio_function_generator.so \
audio_function_generator_wrap.cxx -lstdc++ -lm -lluajit
