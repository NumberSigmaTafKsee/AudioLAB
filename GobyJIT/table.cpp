#include <iostream>
#include <cmath>
#include <vector>

typedef float DspFloatType;

struct BipolarLookupTable
{
    std::vector<DspFloatType> table;
    
    BipolarLookupTable(size_t n) {
        table.resize(n);        
    }
    void sin() {        
        for(int i = -(int)(table.size()/2); i < (int)(table.size()/2); i++)
        {            
            double v = std::sin(2*M_PI*(double)i/((double)table.size()/2.0));            
            table[i+table.size()/2] = v;
        }
    }
    DspFloatType operator()(DspFloatType x) { 
            
        x/=2*M_PI;
        
        DspFloatType r = x * table.size()/2;        
        int i = floor(x);

        DspFloatType f = x - i;

        i += r;        
        if(i >= table.size())  i = table.size()-1;
        
        DspFloatType out;                
        out = table[i] + f*(table[i+1]-table[i]);
        return out;
    }
};

int main()
{
    BipolarLookupTable bt(4096);
    bt.sin();
 
    
    std::cout << bt(1) << ",";
    std::cout << bt(M_PI/2.0) << ",";
    std::cout << bt(M_PI) << ",";
}