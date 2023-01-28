#include <algorithm>
#include <vector>
#include <cmath>
#include <ctime>
#include "SndFile.hpp"

using namespace SoundWave;

void zero(std::vector<float> & x) { memset(x.data(),0x00,sizeof(float)*x.size()); }
void fill(std::vector<float> & x, float v) { for(size_t i = 0; i < x.size(); i++) x[i] = v; }

void circular_convolver(std::vector<float> & a, std::vector<float> &b, std::vector<float> & r) {
    r.resize(a.size()+b.size()-1);    
    zero(r);
    int h_start,x_start,x_end;
    for(int i = 0; i < a.size()+b.size()-1; i++)
    {
        
        x_start = std::max(0,i-(int)a.size()+1);
        x_end   = std::min(i+1,(int)b.size());
        h_start = std::min(i,(int)a.size()-1);
        for(int j = x_start; j < x_end ; j++)
            r[i] += a[h_start--] * b[j];            
    }
}
void linear_convolver(std::vector<float> & a, std::vector<float> &b, std::vector<float> & r) 
{
    r.resize(a.size()+b.size()-1);    
    zero(r);
    for(int i = 0; i < a.size()+b.size()-1; i++)
    {
        r[i] = 0.0f;
            
        int kmin = (i >= b.size() - 1) ? i - (b.size() - 1) : 0;
        int kmax = (i < a.size() - 1) ? i : a.size() - 1;
        for(int j = kmin; j < kmax; j++)
            r[i] = a[i]*(b[i-j]);
    }
}

int main() 
{
    srand(time(NULL));
    SndFileReader a("test.wav");
    SndFileReader b("test.wav");

    std::vector<float> x,y,z;
    x.resize(a.size());
    y.resize(a.size());
    a >> x;    
    b >> y;
    

    linear_convolver(x,y,z);
    
    SndFileWriter w("convolve.wav",a.format(),a.channels(),44100);
    w.write(a.size(),z.data());
    w.close();
}