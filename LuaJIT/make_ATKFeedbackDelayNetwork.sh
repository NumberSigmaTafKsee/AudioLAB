    swig -lua -c++ -Iinclude ATKFeedbackDelayNetwork.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKFeedbackDelayNetwork.so ATKFeedbackDelayNetwork_wrap.cxx -lstdc++ -lm -lluajit    
    