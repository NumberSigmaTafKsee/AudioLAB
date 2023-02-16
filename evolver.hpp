#include <algorithm>
#include <complex>
#include <vector>
#include <climits>
#include <iostream>
#include <random>
#include <chrono>
#include <cstring>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <memory>
#include "StdNoise.hpp"
#include "ga.h"


typedef float floatType;

struct Randomizer
{
    Randomizer() {  }

    static void seed() { srand(time(NULL)); }
    floatType      frand() { return ((floatType)::rand()/(floatType)RAND_MAX); }
    floatType      rand() { return ((floatType)::rand()/(floatType)RAND_MAX); }
    uint64_t    randint(int min, int max) { return round((max-min)*frand())+min; }
    bool        flip(floatType prob) { return frand() < prob; }
    uint64_t    random(int mod) { return ::rand() % mod; }    
};

inline bool coin(double p)
{
    Randomizer noise;
    return noise.rand() < p;
}

inline uint64_t reverse(uint64_t value)
{
	uint64_t r = 0;
	#pragma omp simd
	for(size_t i = 0; i < 64; i++)
	{
		uint64_t bit = (value & (1 << i)) >> i;
		r |= bit >> (64-i);
	}
	return r;
}

inline uint64_t extract_int(unsigned int * p, size_t start, size_t len, int mini=1)
{
	if(len > 64) len = 64;
	uint64_t x = 0;
	#pragma omp simd
	for(size_t i = start; i < (start+len); i++)
	{
		int pos = i/(8*sizeof(unsigned int));
		int bit = i % (8*sizeof(unsigned int));				
		unsigned int q = p[pos];
		x |= (q & (1 <<bit)) >> start;
	}
	return mini*x;
}

inline int8_t extract_int8(unsigned int * p, size_t start, int mini=1)
{
	size_t len = 8;
	int8_t x = 0;
	#pragma omp simd
	for(size_t i = start; i < (start+len); i++)
	{
		int pos = i/(8*sizeof(unsigned int));
		int bit = i % (8*sizeof(unsigned int));		
		unsigned int q = p[pos];
		x |= (q & (1 <<bit));
	}
	return mini*x;
}

inline int16_t extract_int16(unsigned int * p, size_t start, int mini=1)
{
	size_t len = 16;
	int16_t x = 0;
	#pragma omp simd
	for(size_t i = start; i < (start+len); i++)
	{
		int pos = i/(8*sizeof(unsigned int));
		int bit = i % (8*sizeof(unsigned int));		
		unsigned int q = p[pos];
		x |= (q & (1 <<bit));
	}
	return mini*x;
}

inline int32_t extract_int32(unsigned int * p, size_t start, int mini=1)
{
	size_t len = 32;
	int32_t x = 0;
	#pragma omp simd
	for(size_t i = start; i < (start+len); i++)
	{
		int pos = i/(8*sizeof(unsigned int));
		int bit = i % (8*sizeof(unsigned int));		
		unsigned int q = p[pos];
		x |= (q & (1 <<bit));
	}
	return mini*x;
}

inline int64_t extract_int64(unsigned int * p, size_t start, int mini=1)
{
	size_t len = 64;
	int64_t x = 0;
	#pragma omp simd
	for(size_t i = start; i < (start+len); i++)
	{
		int pos = i/(8*sizeof(unsigned int));
		int bit = i % (8*sizeof(unsigned int));		
		unsigned int q = p[pos];
		x |= (q & (1 <<bit));
	}
	return mini*x;
}

inline uint64_t max_bits(size_t len) {
	uint64_t k =  0;
	#pragma omp simd
	for(size_t i = 0; i < len; i++)
		k |= 1 << i;
	return k;
}

inline double extract_double(unsigned * p, size_t start, size_t len, double min, double max)
{
	uint64_t x = extract_int(p,start,len,-1);	
	double  r = (double)x / (double)max_bits(len);		
	return (max-min)*r + min;
}

struct GAChromosome
{
	struct chromosome * p;

	GAChromosome(struct chromosome * c) : p(c) {}
	
	size_t size() const { return p->len; }
	unsigned int * allele() { return p->allele; }

	unsigned int operator[](size_t i) { return p->allele[i]; }
	unsigned int operator()(size_t i) { return chromosome_get_allele(p,i); }
	
	void set(unsigned int pos) { chromosome_set_allele(p,pos); }
	unsigned int get(unsigned int pos) { return chromosome_get_allele(p,pos); }
	void flip(unsigned int pos) { chromosome_not_allele(p,pos); }
	
	std::string asString() { 
		std::string r;
		char * ps = chromosome_as_string(p);
		r = ps;
		free(ps);
		return r;
	}	

	
	// encode paraemeters into the genetic algorithm
	void encode_parameters(void * data)
	{		
		memcpy(p->allele,data,p->len/8);		
	}

	// decode parameters from the genetic algorithm
	template<typename T>
	T* decode_parameters() {
		T * r = (T*)p->allele;
		return r;
	}
	
};

struct GAIndividual
{
	struct individual * p;
	GAIndividual(struct individual * i) : p(i) {}
	GAChromosome get_chromosome() { return GAChromosome(p->chrom); }
	double get_fitness() { return p->fitness; }

	void random() { individual_random(p); }
	int  compare(GAIndividual & i) { return individual_compare(&p,&i.p); }

	bool operator < (GAIndividual& i) { return compare(i) < 0; }
	bool operator > (GAIndividual& i) { return compare(i) > 0; }
	bool operator == (GAIndividual& i) { return compare(i) == 0; }
	bool operator <= (GAIndividual& i) { return compare(i) < 0 || compare(i) == 0; }
	bool operator >= (GAIndividual& i) { return compare(i) > 0 || compare(i) == 0; }

	void print() { individual_print(p,stdout); }
	GAChromosome chromosome() { return GAChromosome(individual_get_chromosome(p)); }
	void setFitness(double v) {  p->fitness = v; }
	double getFitness() const { return p->fitness; }
};

struct GAFittest 
{
	struct fittest * p;
	GAFittest(struct fittest * c) : p(c) {}
	GAIndividual individual() { return GAIndividual(p->i); }
	int generation() const { return p->generation; }
};

struct GAPopulation
{
	struct population * p;

	GAPopulation(struct population * c) : p(c) {}
	size_t size() const { return p->len; }
	GAIndividual operator[](size_t i) { return GAIndividual(p->pop[i]); }
	GAIndividual fittest() const { return GAIndividual(p->fittest); }
	struct fitness_stats stats() const { return p->stats; }
	size_t mutations() const { return p->mutations; }
	size_t crossovers() const { return p->crossovers; }

	void compute_fitness_stats() { population_compute_fitness_stats(p); }
	GAIndividual getFittest() { return GAIndividual(population_get_fittest(p)); }
	void print() { population_print(p,stdout); }
	
};


struct GARandom
{
	GARandom(unsigned long seed) {
		random_seed(seed);
	}
	double random() { return random_random(); }
	int    flip(float p) { return random_flip(p); }
	int    range(int low, int high) { 	return random_range(low,high); }
};

struct GA
{
	struct ga * p;

	GA(const char * filename)
	{
		load(filename);
	}
	GA(unsigned int max_gen,
			size_t chrom_len,
			size_t initial, size_t normal,
			float pcrossover, float pmutation,
			enum ga_selection_strategies selection_strategy,
			enum ga_crossover_strategies crossover_strategy,
			objective_fn obj_fn) {
		
		p = new_ga(max_gen,
			chrom_len,
			initial,normal,
			pcrossover,pmutation,
			selection_strategy,
			crossover_strategy,
			obj_fn);
	}
	
	~GA() {
		if(p) delete_ga(p);
	}

	GAIndividual operator[](size_t i) {
		return pop()[i];
	}

	unsigned current() const { return p->current; }
	unsigned maxgen()  const { return p->max_gen; }
	size_t   chromlen() const { return p->chrom_len; }
	size_t   initial() const { return p->initial; }
	size_t   normal()  const { return p->normal; }
	float    pcrossover() const { return p->pcrossover; }
	float    pmutation() const  { return p->pmutation; }

	GAPopulation pop() { return GAPopulation(p->cur_pop); }
	GAPopulation old() { return GAPopulation(p->old_pop); }
	GAFittest    best() { return GAFittest(&p->best); }

	void set_report_strategy(enum ga_report_strategies report_strategy) {
		ga_set_report_strategy(p,report_strategy);
    }
	
	void evolve(unsigned maxgen) {
		ga_evolve(p,maxgen);
	}
	void reset() {
		p->current = 0;
	}	
	void step() {
		ga_step(p);
	}	
	GAFittest get_best_ever() { return GAFittest(ga_get_best_ever(p)); }

	void save(const char * filename)
	{

	}
	void load(const char * filename)
	{

	}
};

std::vector<floatType> examples;
std::vector<floatType> training;
