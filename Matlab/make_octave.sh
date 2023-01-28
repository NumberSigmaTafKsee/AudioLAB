PKG_CONFIG_PATH=~/pkgconfig gcc `pkg-config --cflags octave octinterp`  -o octave octave.cpp -lstdc++ -lm  `pkg-config --libs octave octinterp`
