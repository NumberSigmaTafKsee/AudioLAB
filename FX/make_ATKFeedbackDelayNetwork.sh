    swig -lua -c++ -I../include -I../include ATKFeedbackDelayNetwork.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKFeedbackDelayNetwork.so ATKFeedbackDelayNetwork_wrap.cxx -lstdc++ -lm -lluajit 
    