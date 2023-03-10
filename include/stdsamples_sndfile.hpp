#pragma once 

extern "C" {
#include <sndfile.h>
#include <stdint.h>
}

#include <vector>
#include <cassert>
#include <cstring>


namespace AudioDSP {
	template<typename T>
	struct SndFileReader
	{
		SNDFILE * sndfile;
		SF_INFO   info;

		SndFileReader(const char * filename) {
			sndfile = sf_open(filename,SFM_READ,&info);
			assert(sndfile != NULL);
		}
		~SndFileReader() {
			if(sndfile) sf_close(sndfile);
		}    

		int64_t size() const {
			return  frames() * channels();
		}
		
		virtual SndFileReader<T>& operator >> (std::vector<T> & s) = 0;

		virtual int64_t read(size_t n, std::vector<T> & samples) = 0;
		virtual int64_t read(size_t n, T * samples) = 0;

		virtual int64_t read_frames(size_t n, std::vector<T> & samples) = 0;
		virtual int64_t read_frames(size_t n, T * samples) = 0;

		int64_t seek(int64_t frames, int whence = SEEK_SET) {
			return sf_seek(sndfile,frames,whence);
		}

		const char* get_string(int str_type) {
			return sf_get_string(sndfile,str_type);
		}
		int command(int cmd, void * data, int datasize) {
			return sf_command(sndfile,cmd,data,datasize);
		}

		int64_t frames() const { return info.frames; }
		int samplerate() const { return info.samplerate; }
		int channels() const { return info.channels; }
		int format() const { return info.format;}
		int sections() const { return info.sections; }
		int seekable() const { return info.seekable; }
		

	};

	struct SndFileReaderFloat : public SndFileReader<float>
	{
		SndFileReaderFloat(const char * filename) 
		: SndFileReader<float>(filename)
		{
		
		}

		using SndFileReader::size;
		using SndFileReader::seek;
		using SndFileReader::get_string;
		using SndFileReader::command;
		using SndFileReader::frames;
		using SndFileReader::samplerate;
		using SndFileReader::channels;
		using SndFileReader::format;
		using SndFileReader::sections;
		using SndFileReader::seekable;

		int64_t read(size_t n, std::vector<float> & samples) {        
			samples.resize(n);
			return sf_read_float(sndfile,samples.data(),n);
		}
		int64_t read(size_t n, float * samples) {                    
			return sf_read_float(sndfile,samples,n);
		}

		SndFileReader<float>& operator >> (std::vector<float> & s) {
			sf_read_float(sndfile,s.data(),s.size());
			return *this;
		}

		int64_t read_frames(size_t n, std::vector<float> & samples) {        
			samples.resize(n);
			return sf_readf_float(sndfile,samples.data(),n);
		}
		int64_t read_frames(size_t n, float * samples) {                    
			return sf_readf_float(sndfile,samples,n);
		}

	};

	struct SndFileReaderDouble : public SndFileReader<double>
	{
		SndFileReaderDouble(const char * filename) 
		: SndFileReader<double>(filename)
		{
		
		}

		using SndFileReader::size;
		using SndFileReader::seek;
		using SndFileReader::get_string;
		using SndFileReader::command;
		using SndFileReader::frames;
		using SndFileReader::samplerate;
		using SndFileReader::channels;
		using SndFileReader::format;
		using SndFileReader::sections;
		using SndFileReader::seekable;

		SndFileReader<double>& operator >> (std::vector<double> & s) {
			sf_read_double(sndfile,s.data(),s.size());
			return *this;
		}
		int64_t read(size_t n, std::vector<double> & samples) {        
			samples.resize(n);
			return sf_read_double(sndfile,samples.data(),n);
		}
		int64_t read(size_t n, double * samples) {                    
			return sf_read_double(sndfile,samples,n);
		}

		int64_t read_frames(size_t n, std::vector<double> & samples) {        
			samples.resize(n);
			return sf_readf_double(sndfile,samples.data(),n);
		}
		int64_t read_frames(size_t n, double * samples) {                    
			return sf_readf_double(sndfile,samples,n);
		}

	};

	struct SndFileReaderInt : public SndFileReader<int>
	{
		SndFileReaderInt(const char * filename) 
		  : SndFileReader<int>(filename)
		{
		
		}

		using SndFileReader::size;
		using SndFileReader::seek;
		using SndFileReader::get_string;
		using SndFileReader::command;
		using SndFileReader::frames;
		using SndFileReader::samplerate;
		using SndFileReader::channels;
		using SndFileReader::format;
		using SndFileReader::sections;
		using SndFileReader::seekable;

		SndFileReader<int>& operator >> (std::vector<int> & s) {
			sf_read_int(sndfile,s.data(),s.size());
			return *this;
		}
		int64_t read(size_t n, std::vector<int> & samples) {        
			samples.resize(n);
			return sf_read_int(sndfile,samples.data(),n);
		}
		int64_t read(size_t n, int * samples) {                    
			return sf_read_int(sndfile,samples,n);
		}

		int64_t read_frames(size_t n, std::vector<int> & samples) {        
			samples.resize(n);
			return sf_readf_int(sndfile,samples.data(),n);
		}
		int64_t read_frames(size_t n, int * samples) {                    
			return sf_readf_int(sndfile,samples,n);
		}

	};

	struct SndFileReaderShort : public SndFileReader<short>
	{
		
		SndFileReaderShort(const char * filename) 
		: SndFileReader<short>(filename)
		{
		
		}

		using SndFileReader::size;
		using SndFileReader::seek;
		using SndFileReader::get_string;
		using SndFileReader::command;
		using SndFileReader::frames;
		using SndFileReader::samplerate;
		using SndFileReader::channels;
		using SndFileReader::format;
		using SndFileReader::sections;
		using SndFileReader::seekable;

		SndFileReader<short>& operator >> (std::vector<short> & s) {
			sf_read_short(sndfile,s.data(),s.size());
			return *this;
		}
		int64_t read(size_t n, std::vector<short> & samples) {        
			samples.resize(n);
			return sf_read_short(sndfile,samples.data(),n);
		}
		int64_t read(size_t n, short * samples) {                    
			return sf_read_short(sndfile,samples,n);
		}
		int64_t read_frames(size_t n, std::vector<short> & samples) {        
			samples.resize(n);
			return sf_readf_short(sndfile,samples.data(),n);
		}
		int64_t read_frames(size_t n, short * samples) {                    
			return sf_readf_short(sndfile,samples,n);
		}
	};



	template<typename T>
	struct SndFileWriter
	{
		SNDFILE * sndfile;
		SF_INFO   info;

		SndFileWriter(const char * filename, int format, int channels, int sample_rate) {
			memset(&info,0x00,sizeof(info));
			info.channels = channels;
			info.format = format;
			info.samplerate = sample_rate;
			sndfile = sf_open(filename,SFM_WRITE,&info);
		}
		~SndFileWriter() {
			if(sndfile) sf_close(sndfile);
		}
		void close() {
			if(sndfile) {
				sf_close(sndfile);
				sndfile = NULL;
			}
		}

		
		virtual SndFileWriter<T>& operator << (std::vector<T> & s) = 0;

		virtual void write(std::vector<T> & v) = 0;
		virtual void write(size_t n,T * v) = 0;

		virtual void write_frames(std::vector<T> & v) = 0;
		virtual void write_frames(size_t n, T * v) = 0;

		
		int64_t seek(int64_t frames, int whence = SEEK_SET) {
			return sf_seek(sndfile,frames,whence);
		}
		int set_string(int str_type, const char * str) {
			return sf_set_string(sndfile,str_type,str);
		}
		void sync() {
			sf_write_sync(sndfile);
		}
		
		int command(int cmd, void * data, int datasize) {
			return sf_command(sndfile,cmd,data,datasize);
		}

		int samplerate() const { return info.samplerate; }
		int channels() const { return info.channels; }
		int format() const { return info.format;}
		int sections() const { return info.sections; }
		int seekable() const { return info.seekable; }

	};


	struct SndFileWriterFloat : public SndFileWriter<float>
	{
		SndFileWriterFloat(const char * filename, int format, int channels, int sample_rate)
		: SndFileWriter<float>(filename,format,channels,sample_rate)
		{

		}

		SndFileWriter<float>& operator << (std::vector<float> & s) {
			sf_write_float(sndfile,s.data(),s.size());
			return *this;
		}
		
		
		using SndFileWriter::seek;
		using SndFileWriter::sync;
		using SndFileWriter::set_string;
		using SndFileWriter::command;    
		using SndFileWriter::samplerate;
		using SndFileWriter::channels;
		using SndFileWriter::format;
		using SndFileWriter::sections;
		using SndFileWriter::seekable;

		void write(std::vector<float> & v) {
			sf_write_float(sndfile,v.data(),v.size());
		}
		void write(size_t n,float * v) {
			sf_write_float(sndfile,v,n);
		}
		void write_frames(std::vector<float> & v) {
			sf_writef_float(sndfile,v.data(),v.size());
		}
		void write_frames(size_t n, float * v) {
			sf_writef_float(sndfile,v,n);
		}

	};

	struct SndFileWriterDouble : public SndFileWriter<double>
	{

		SndFileWriterDouble(const char * filename, int format, int channels, int sample_rate)
		: SndFileWriter<double>(filename,format,channels,sample_rate)
		{

		}
		
		SndFileWriter<double>& operator << (std::vector<double> & s) {
			sf_write_double(sndfile,s.data(),s.size());
			return *this;
		}

		
		using SndFileWriter::seek;
		using SndFileWriter::sync;
		using SndFileWriter::set_string;
		using SndFileWriter::command;    
		using SndFileWriter::samplerate;
		using SndFileWriter::channels;
		using SndFileWriter::format;
		using SndFileWriter::sections;
		using SndFileWriter::seekable;

		void write(std::vector<double> & v) {
			sf_write_double(sndfile,v.data(),v.size());
		}
		void write(size_t n,double * v) {
			sf_write_double(sndfile,v,n);
		}
		void write_frames(std::vector<double> & v) {
			sf_writef_double(sndfile,v.data(),v.size());
		}
		void write_frames(size_t n,double * v) {
			sf_writef_double(sndfile,v,n);
		}

	};

	struct SndFileWriterShort : public SndFileWriter<short>
	{
		SndFileWriterShort(const char * filename, int format, int channels, int sample_rate)
		: SndFileWriter<short>(filename,format,channels,sample_rate)
		{

		}
		
		
		using SndFileWriter::seek;
		using SndFileWriter::sync;
		using SndFileWriter::set_string;
		using SndFileWriter::command;    
		using SndFileWriter::samplerate;
		using SndFileWriter::channels;
		using SndFileWriter::format;
		using SndFileWriter::sections;
		using SndFileWriter::seekable;

		SndFileWriter<short>& operator << (std::vector<short> & s) {
			sf_write_short(sndfile,s.data(),s.size());
			return *this;
		}
		void write(std::vector<short> & v) {
			sf_write_short(sndfile,v.data(),v.size());
		}
		void write(size_t n,short * v) {
			sf_write_short(sndfile,v,n);
		}
		void write_frames(std::vector<short> & v) {
			sf_writef_short(sndfile,v.data(),v.size());
		}
		void write_frames(size_t n,short * v) {
			sf_writef_short(sndfile,v,n);
		}

	};


	struct SndFileWriterInt : public SndFileWriter<int>
	{
		SndFileWriterInt(const char * filename, int format, int channels, int sample_rate)
		: SndFileWriter<int>(filename,format,channels,sample_rate)
		{

		}
		
		
		using SndFileWriter::seek;
		using SndFileWriter::sync;
		using SndFileWriter::set_string;
		using SndFileWriter::command;    
		using SndFileWriter::samplerate;
		using SndFileWriter::channels;
		using SndFileWriter::format;
		using SndFileWriter::sections;
		using SndFileWriter::seekable;

		SndFileWriter<int>& operator << (std::vector<int> & s) {
			sf_write_int(sndfile,s.data(),s.size());
			return *this;
		}
		void write(std::vector<int> & v) {
			sf_write_int(sndfile,v.data(),v.size());
		}
		void write(size_t n,int * v) {
			sf_write_int(sndfile,v,n);
		}
		void write_frames(std::vector<int> & v) {
			sf_writef_int(sndfile,v.data(),v.size());
		}
		void write_frames(size_t n,int * v) {
			sf_writef_int(sndfile,v,n);
		}

	};
}
