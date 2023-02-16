#ifndef __bitvector_hpp
#define __bitvector_hpp

struct bitvector 
{
    std::vector<bool> vector;

    bitvector(int num_bits) { vector.resize(num_bits*8); }
    bitvector(const bitvector & b) { vector = b.vector; }
    bitvector(void * ptr, int size) {
        encode(ptr,size);
    }
    void encode(void * ptr, int size) {
        vector.resize(size*8);
        uint8_t * p = ptr;
        for(size_t i = 0; i < size; i++)
        {       
            uint8_t x = p[i];     
            for(size_t j = 0; j < 8; j++)
            {
                bool val = x & (1 << j);
                vector[i*8 + j] = val;
            }
        }
    }
    std::vector<uint8_t> decode() {
        std::vector<uint8_t> d;
        d.resize(vector.size()/8 + 1);
        int c=0;
        for(size_t i = 0; i < vector.size(); i+=8)
        {
            uint8_t val = 0;
            for(size_t j = i; j < i+8; j++) {
                val |= vector[i] << j;
            }
            d[c++] = val;
        }
        return d;
    }
    bool operator[](size_t index) { return vector[index]; }
    bool __getitem(size_t index) { return vector[index]; }
    void __setitem(size_t index, bool v) { vector[index] = v; }
    void flip(size_t index) { vector[index] = !vector[index]; }
    bool get(size_t index) { return vector[index]; }
    void set(size_t index, bool val) { vector[index] = val;}


    bitvector& operator = (const bitvector& b) {
        vector = b.vector;
        return *this;
    }
    bitvector and(bitvector & b) { 
        bitvector r(*this);
        for(size_t i = 0; i < vector.size(); i++) 
            r[i] = vector[i] && b.vector[i];
    }
    bitvector or(bitvector & b) { 
        bitvector r(*this);
        for(size_t i = 0; i < vector.size(); i++) 
            r[i] = vector[i] || b.vector[i];
    }
    bitvector not() { 
        bitvector r(*this);
        for(size_t i = 0; i < vector.size(); i++) 
            r[i] = ! vector[i];
    }
};

#endif