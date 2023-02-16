%module file
%{
#include <cstdio>
#include <string>
#include <iostream>
#include <vector>

%}

%include "stdint.i"
%include "std_string.i"

%inline %{
     enum SeekType
 {
     BEGIN = SEEK_SET,
     CURRENT= SEEK_CUR,
     END = SEEK_END,
 };

 class StdFile
 {
    protected:

        FILE * file;

    public:

    
        StdFile(const char * filename, const char * mode) {
            file = fopen(filename,mode);
            assert(file != NULL);
        }
        ~StdFile() { if(file) fclose(file); }

        void open(const char * filename, const char * mode) {
            if(file) fclose(file);
            file = fopen(filename,mode);
            assert(file != NULL);
        }
        void close() { if(file) { fclose(file); file = NULL;} }

        void seek(uint64_t pos, SeekType w = BEGIN) {
            fseek(file,pos,w);
        }

    };

    class BinaryFileWriter : public StdFile
    {
    public:

        BinaryFileWriter(const char * file, const char * mode = "wb") : StdFile(file,mode) {        
        }
        
        void write_i8(int8_t i) { fwrite(&i,sizeof(i),1,file); }
        void write_i16(int16_t i) { fwrite(&i,sizeof(i),1,file); }
        void write_i32(int32_t i) { fwrite(&i,sizeof(i),1,file); }
        void write_i64(int64_t i) { fwrite(&i,sizeof(i),1,file); }        
        void write_u8(uint8_t i) { fwrite(&i,sizeof(i),1,file); }
        void write_u16(uint16_t i) { fwrite(&i,sizeof(i),1,file); }
        void write_u32(uint32_t i) { fwrite(&i,sizeof(i),1,file); }
        void write_u64(uint64_t i) { fwrite(&i,sizeof(i),1,file); }
        void write_str(const char * s) {
            size_t len = strlen(s);
            write(len);
            fwrite(s,1,strlen(s),file);
        }
        void write_string(const std::string& s) {
            write(s.data());
        }

        template<typename T>
        void write(T i) { fwrite(&i,sizeof(T),1,file); }

        template<typename T>
        void write(const std::vector<T> & v)
        {
            size_t len = v.size();
            write(len);
            fwrite(v.data(),sizeof(T),v.size(),file);
        }
                
        
        void write_i8vec(const std::vector<int8_t> & v) { write<int8_t>(v); }
        void write_u8vec(const std::vector<uint8_t> & v) { write<uint8_t>(v); }
        void write_i16vec(const std::vector<int16_t> & v) { write<int16_t>(v); }
        void write_u16vec(const std::vector<uint16_t> & v) { write<uint16_t>(v); }
        void write_i32vec(const std::vector<int32_t> & v) { write<int32_t>(v); }
        void write_u32vec(const std::vector<uint32_t> & v) { write<uint32_t>(v); }
        void write_i64vec(const std::vector<int64_t> & v) { write<int64_t>(v); }
        void write_u64vec(const std::vector<uint64_t> & v) { write<uint64_t>(v); }
        void write_fvec(const std::vector<float> & v) { write<float>(v); }
        void write_dvec(const std::vector<double> & v) { write<double>(v); }        
    };

    class BinaryFileReader : public StdFile
    {
    public:

        BinaryFileReader(const char * file, const char * mode = "rb") : StdFile(file,mode) {        
        }

        
        size_t read_i8(int8_t &i)  { return fread(&i,sizeof(i),1,file); }
        size_t read_i16(int16_t &i) { return fread(&i,sizeof(i),1,file); }
        size_t read_i32(int32_t &i) { return fread(&i,sizeof(i),1,file); }
        size_t read_i64(int64_t &i) { return fread(&i,sizeof(i),1,file); }        
        size_t read_u8(uint8_t &i) { return fread(&i,sizeof(i),1,file); }
        size_t read_u16(uint16_t &i) { return fread(&i,sizeof(i),1,file); }
        size_t read_u32(uint32_t &i) { return fread(&i,sizeof(i),1,file); }
        size_t read_u64(uint64_t &i) { return fread(&i,sizeof(i),1,file); }
                        
        
        template<typename T>
        size_t read(T i) { return fread(&i,sizeof(T),1,file); }

        std::string read_string() {
            size_t n=0;
            std::string s;
            read(n);
            if(n > 0) {
                s.resize(n+1);
                memset((void*)s.data(),0,n+1);
                fread((void*)s.data(),1,n,file);
            }
            return s;
        }
        
        template<typename T>
        size_t read(std::vector<T> & v)
        {
            size_t n = 0;
            read(n);            
            if(n > 0) {
                v.resize(n);
                n = fread(v.data(),sizeof(T),v.size(),file);
            }
            return n;
        }

        size_t read_i8vec(const std::vector<int8_t> & v) { return read(v); }
        size_t read_u8vec(const std::vector<uint8_t> & v) { return read(v); }
        size_t read_i16vec(const std::vector<int16_t> & v) { return read(v); }
        size_t read_u16vec(const std::vector<uint16_t> & v) { return read(v); }
        size_t read_i32vec(const std::vector<int32_t> & v) { return read(v); }
        size_t read_u32vec(const std::vector<uint32_t> & v) { return read(v); }
        size_t read_i64vec(const std::vector<int64_t> & v) { return read(v); }
        size_t read_u64vec(const std::vector<uint64_t> & v) { return read(v); }
        size_t read_fvec(const std::vector<float> & v) { return read(v); }
        size_t read_dvec(const std::vector<double> & v) { return read(v); }                
    };
%}