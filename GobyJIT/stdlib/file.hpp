#ifndef __FILE_HPP
#define __FILE_HPP

struct FileWriter
{
    FILE * f;

    FileWriter(const char * filename, bool append=false) {
        std::string p = "wb";        
        if(append) p += "+";
        f = fopen(filename,p.c_str());
        assert(f != NULL);
    }
    ~FileWriter() {
        if(f != NULL) fclose(f);
    }

    void write_i8(int8_t val) {
        fwrite(&val,sizeof(int8_t),1,f);
    }
    void write_i16(int16_t val) {
        fwrite(&val,sizeof(int16_t),1,f);
    }
    void write_i32(int32_t val) {
        fwrite(&val,sizeof(int32_t),1,f);
    }
    void write_i64(int32_t val) {
        fwrite(&val,sizeof(int64_t),1,f);
    }
    void write_ui8(uint8_t val) {
        fwrite(&val,sizeof(uint8_t),1,f);
    }
    void write_ui16(uint16_t val) {
        fwrite(&val,sizeof(uint16_t),1,f);
    }
    void write_ui32(uint32_t val) {
        fwrite(&val,sizeof(uint32_t),1,f);
    }
    void write_ui64(uint64_t val) {
        fwrite(&val,sizeof(uint64_t),1,f);
    }
    void write_float(float val) {
        fwrite(&val,sizeof(float),1,f);
    }
    void write_double(uint64_t val) {
        fwrite(&val,sizeof(double),1,f);
    }
    void write_string(const char * s) {        
        write_u32(strlen(s));
        fwrite(s,sizeof(char),strlen(s),f);
    }
    void write_size(size_t val) {
        fwrite(&val,sizeof(size_t),1,f);
    }    
    void write_data(void * ptr, size_t size) {
        write_size(size);
        fwrite(ptr,1,size,f);
    }
    void write_vi8(const std::vector<int8_t> & v) {
        write_size(v.size());
        fwrite(v.data(),sizeof(int8_t),v.size(),1);
    }    
    void write_vi16(const std::vector<int16_t> & v) {
        write_size(v.size());
        fwrite(v.data(),sizeof(int16_t),v.size(),1);
    }    
    void write_vi32(const std::vector<int32_t> & v) {
        write_size(v.size());
        fwrite(v.data(),sizeof(int32_t),v.size(),1);
    }    
    void write_vi64(const std::vector<int64_t> & v) {
        write_size(v.size());
        fwrite(v.data(),sizeof(int64_t),v.size(),1);
    }

    void write_vui8(const std::vector<uint8_t> & v) {
        write_size(v.size());
        fwrite(v.data(),sizeof(uint8_t),v.size(),1);
    }    
    void write_vui16(const std::vector<uint16_t> & v) {
        write_size(v.size());
        fwrite(v.data(),sizeof(uint16_t),v.size(),1);
    }    
    void write_vui32(const std::vector<uint32_t> & v) {
        write_size(v.size());
        fwrite(v.data(),sizeof(uint32_t),v.size(),1);
    }    
    void write_vui64(const std::vector<uint64_t> & v) {
        write_size(v.size());
        fwrite(v.data(),sizeof(uint64_t),v.size(),1);
    }

    void write_vfloat(const std::vector<float> & v) {
        write_size(v.size());
        fwrite(v.data(),sizeof(float),v.size(),1);
    }    
    void write_vdouble(const std::vector<double> & v) {
        write_size(v.size());
        fwrite(v.data(),sizeof(double),v.size(),1);        
    }    
    
    template<typename T>
    void write_vector(std::vector<T> & v) {
        write_size(v.size());
        fwrite(v.data(),sizeof(T),v.size(),f);
    }
};


struct FileReader 
{
    FILE * f;

    FileReader(const char * filename) { f = fopen(f,"rb"); assert(f != NULL); }
    ~FileReader() { if(f) fclose(f); }

    template<typename T>
    void read(T& val) {
        fread(&val,sizeof(T),1,f);
    }
    template<typename T>
    void read_vector(std::vector<T> & v) {
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        v.resize(i);
        fread(v.data(),sizeof(T),i,f);
    }
    
    int8_t read_i8() { int8_t r; fread(&r,sizeof(int8_t),1,f); return r; }
    uint8_t read_ui8() { uint8_t r; fread(&r,sizeof(uint8_t),1,f); return r; }
    int16_t read_i16() { int16_t r; fread(&r,sizeof(int16_t),1,f); return r; }
    uint16_t read_ui16() { uint16_t r; fread(&r,sizeof(uint16_t),1,f); return r; }
    int32_t read_i32() { int32_t r; fread(&r,sizeof(int32_t),1,f); return r; }
    uint32_t read_ui32() { uint32_t r; fread(&r,sizeof(uint32_t),1,f); return r; }
    int64_t read_i64() { int64_t r; fread(&r,sizeof(int64_t),1,f); return r; }
    uint64_t read_ui64() { uint64_t r; fread(&r,sizeof(uint64_t),1,f); return r; }
    float read_float() { float r; fread(&r,sizeof(float),1,f); return r; }
    double read_double() { double r; fread(&r,sizeof(double),1,f); return r; }
    size_t read_size() { size_t s; fread(&s,sizeof(size_t),1,f); return s; }

    std::string read_i8() { 
        std::string r; 
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        char buffer[i+1];
        memset(buffer,0x00,i);
        fread(buffer,sizeof(char),i,f); return r; }
        r = buffer;
        return r;
    }    

    std::vector<int8_t> read_vi8() { 
        std::vector<int8_t> r;
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        r.resize(i);
        fread(r.data(),sizeof(int8_t),i,f);
        return r;
    }
    std::vector<uint8_t> read_vui8() { 
        std::vector<uint8_t> r;
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        r.resize(i);
        fread(r.data(),sizeof(uint8_t),i,f);
        return r;
    }
    std::vector<int16_t> read_vi16() { 
        std::vector<int8_t> r;
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        r.resize(i);
        fread(r.data(),sizeof(int16_t),i,f);
        return r;
    }
    std::vector<uint16_t> read_vui16() { 
        std::vector<uint16_t> r;
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        r.resize(i);
        fread(r.data(),sizeof(uint16_t),i,f);
        return r;
    }
    std::vector<int32_t> read_vi8() { 
        std::vector<int32_t> r;
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        r.resize(i);
        fread(r.data(),sizeof(int32_t),i,f);
        return r;
    }
    std::vector<uint32_t> read_vui32() { 
        std::vector<uint32_t> r;
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        r.resize(i);
        fread(r.data(),sizeof(uint32_t),i,f);
        return r;
    }
    std::vector<int64_t> read_vi8() { 
        std::vector<int64_t> r;
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        r.resize(i);
        fread(r.data(),sizeof(int64_t),i,f);
        return r;
    }
    std::vector<uint64_t> read_vui64() { 
        std::vector<uint64_t> r;
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        r.resize(i);
        fread(r.data(),sizeof(uint64_t),i,f);
        return r;        
    }
    std::vector<float> read_vfloat() { 
        std::vector<float> r;
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        r.resize(i);
        fread(r.data(),sizeof(float),i,f);
        return r;
    }
    std::vector<uint8_t> read_vdouble() { 
        std::vector<double> r;
        size_t i;
        fread(&i,sizeof(size_t),1,f);
        r.resize(i);
        fread(r.data(),sizeof(double),i,f);
        return r;
    }

    void* read_data() {
        size_t i;
        fread(&i,sizeof(size_t),1);
        void * ptr = calloc(i,sizeof(uint8_t));
        fread(ptr,1,i,f);
        return ptr;
    }
};

#endif
