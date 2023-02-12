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
#include "Plot.hpp"
Default noise;

Plot<float> plot;

bool coin(double p)
{
    return noise.rand() < p;
}

uint64_t reverse(uint64_t value)
{
	uint64_t r = 0;
	for(size_t i = 0; i < 64; i++)
	{
		uint64_t bit = (value & (1 << i)) >> i;
		r |= bit >> (64-i);
	}
	return r;
}

uint64_t extract_int(unsigned int * p, size_t start, size_t len, int mini=1)
{
	if(len > 64) len = 64;
	uint64_t x = 0;
	for(size_t i = start; i < (start+len); i++)
	{
		int pos = i/(8*sizeof(unsigned int));
		int bit = i % (8*sizeof(unsigned int));				
		unsigned int q = p[pos];
		x |= (q & (1 <<bit)) >> start;
	}
	return mini*x;
}

int8_t extract_int8(unsigned int * p, size_t start, int mini=1)
{
	size_t len = 8;
	int8_t x = 0;
	for(size_t i = start; i < (start+len); i++)
	{
		int pos = i/(8*sizeof(unsigned int));
		int bit = i % (8*sizeof(unsigned int));		
		unsigned int q = p[pos];
		x |= (q & (1 <<bit));
	}
	return mini*x;
}
int16_t extract_int16(unsigned int * p, size_t start, int mini=1)
{
	size_t len = 16;
	int16_t x = 0;
	for(size_t i = start; i < (start+len); i++)
	{
		int pos = i/(8*sizeof(unsigned int));
		int bit = i % (8*sizeof(unsigned int));		
		unsigned int q = p[pos];
		x |= (q & (1 <<bit));
	}
	return mini*x;
}
int32_t extract_int32(unsigned int * p, size_t start, int mini=1)
{
	size_t len = 32;
	int32_t x = 0;
	for(size_t i = start; i < (start+len); i++)
	{
		int pos = i/(8*sizeof(unsigned int));
		int bit = i % (8*sizeof(unsigned int));		
		unsigned int q = p[pos];
		x |= (q & (1 <<bit));
	}
	return mini*x;
}
int64_t extract_int64(unsigned int * p, size_t start, int mini=1)
{
	size_t len = 64;
	int64_t x = 0;
	for(size_t i = start; i < (start+len); i++)
	{
		int pos = i/(8*sizeof(unsigned int));
		int bit = i % (8*sizeof(unsigned int));		
		unsigned int q = p[pos];
		x |= (q & (1 <<bit));
	}
	return mini*x;
}



uint64_t max_bits(size_t len) {
	uint64_t k =  0;
	for(size_t i = 0; i < len; i++)
		k |= 1 << i;
	return k;
}
double extract_double(unsigned * p, size_t start, size_t len, double min, double max)
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
	/*
	void step() {
		ga_step(p);
	}
	*/
	GAFittest get_best_ever() { return GAFittest(ga_get_best_ever(p)); }

	void save(const char * filename)
	{

	}
	void load(const char * filename)
	{

	}
};

typedef float floatType;
std::vector<floatType> examples;
std::vector<floatType> training;

floatType chebyshev(floatType x, unsigned int A[], int order)
{
	// To = 1
	// T1 = x
	// Tn = 2.x.Tn-1 - Tn-2
	// out = sum(Ai*Ti(x)) , i C {1,..,order}
	floatType Tn_2 = 1.0f;
	floatType gibs = (floatType)A[0]/(floatType)UINT_MAX;
	floatType Tn_1 = x*gibs;
	floatType Tn;
	floatType out = Tn_1;

	for(int n=2;n<=order;n++)
	{
		Tn    	=   2.0*x*Tn_1 - Tn_2;
		out    +=   Tn;
		Tn_2 	=   Tn_1;
		Tn_1 	=   (floatType)A[n-1]/(floatType)UINT_MAX;
	}
	return out;
}


///////////////////////////////////////////////////////////////
// Test area
///////////////////////////////////////////////////////////////
static void test_ga(void);
static void objective(struct individual *i);
static double chromosome_to_double(struct chromosome *c);
int check_bits(struct chromosome * c, const char * msg);

#define CHEBYS 8
struct Cheby
{
	unsigned int A[CHEBYS];
};

int compf(const void * a, const void * b) {
	return *(unsigned int*)a - *(unsigned int*)b;
}
void objective(struct individual *i)
{
	double value=0;
	struct chromosome *c;	
	assert(i);	
	c = individual_get_chromosome(i);
	assert(c);		
	Cheby *p = (Cheby*)c->allele;
	std::vector<float> out(128);		
	//qsort(p->A,CHEBYS,sizeof(unsigned int),compf);
	for(size_t i = 1; i < CHEBYS; i++)
		if( ((floatType)p->A[i]/(floatType)UINT_MAX) < sqrt((floatType)p->A[i-1]/(floatType)UINT_MAX)) value += i * 10;
 	for(size_t i = 0; i < 128; i++) {
		auto bobo = chebyshev(examples[i],p->A,CHEBYS) - chebyshev(0,p->A,CHEBYS);
		auto x = std::clamp(bobo,-1.0f,1.0f);
		out[i] = x;
	}	
	
	for(size_t i = 0; i < 128; i++) {
		auto x = fabs(out[i] - training[i]) < 1e-2;
		value += x;
	}	
 	individual_set_fitness(i, value);
}


inline float udo1(float x, float g = 1.0)
{
    float ax = fabs(x*g);
    if(ax == 0) return 0;
    return std::clamp(((x/ax)*(1-exp(g*(x*x)/ax))),-1.0,1.0);
}

void test_cheby()
{
	
    for(int i = 0; i < 128; i++) {
		auto sample = std::sin(2*M_PI*(double)i/128.0);
        examples.push_back(sample);
        training.push_back(std::clamp(2*sample,-1.0,1.0)/2.0);
    }
	plot.setstyle("lines");
    plot.plot_x(examples.data(),examples.size(),"examples");
	plot.plot_x(training.data(),training.size(),"training");
	GA ga(500, sizeof(Cheby)*8, 50, 50, 0.9, 1e-3,GA_S_ROULETTE_WHEEL, GA_X_SINGLE_POINT, objective);
	ga.set_report_strategy(GA_R_HUMAN_READABLE);
	ga.evolve(100);
	Cheby * p = (Cheby*)(ga.p->best.i->chrom->allele);
	std::vector<float> out(128);
	for(size_t i = 0; i < CHEBYS; i++) std::cout << (floatType)p->A[i]/(floatType)UINT_MAX << ",";
	std::cout << std::endl;
	
	for(size_t i = 0; i < 128; i++)
	{
		auto bobo = chebyshev(examples[i],p->A,CHEBYS) - chebyshev(0,p->A,CHEBYS);
		out[i] = std::clamp(bobo,-1.0f,1.0f)/2.0f;
	}
		
	plot.plot_x(out.data(),out.size(),"output");
	
	sleep(1000);
}


int
main(int argc, char *argv[])
{
	srand(time(NULL));
	test_cheby();
	
	exit(EXIT_SUCCESS);
}
/*
struct Test
{
	double x[2];
};

void test_ga(void)
{
	struct ga *g;

	random_seed(time(NULL));

	g = new_ga(100, 88, 500, 100, 0.995, 1e-6, 
			GA_S_ROULETTE_WHEEL, GA_X_SINGLE_POINT, objective);
	assert(g);

	//ga_set_report_strategy(g, GA_R_GRAPH);
	ga_set_report_strategy(g, GA_R_HUMAN_READABLE);

	ga_evolve(g, 10);
	//int x = ((Test*)g->best.i->chrom->allele)->x;
	//int y = ((Test*)g->best.i->chrom->allele)->y;
	//char * kaka = ((Test*)g->best.i->chrom->allele)->kaka;
	//double james = ((Test*)g->best.i->chrom->allele)->james/(double)ULLONG_MAX;
	//james = 10.*james - 5.0;
	//printf("best %d=%f\n x=%d y=%d\n",g->best.generation,g->best.i->fitness,x,y);
	//printf("best %d=%f\n kaka=%s\n",g->best.generation,g->best.i->fitness,kaka);
	//printf("best %d=%f\n james=%f\n",g->best.generation,g->best.i->fitness,james);
    
	delete_ga(g);
}


void objective(struct individual *i)
{
	double value=0;
	struct chromosome *c;
	static const char * msg = "Hello World";
	assert(i);
	Test * t;

	c = individual_get_chromosome(i);
	assert(c);
	t = (Test*)c->allele;
	//value = t->x + t->y;	
	//for(size_t i = 0; i < 11; i++)
	//	if(t->kaka[i] == msg[i]) value++;
	value = chromosome_to_double(c)/(double)ULLONG_MAX;
	value = 10.0*value - 5.0;
	individual_set_fitness(i, value );
}

int check_bits(struct chromosome * c, const char * msg)  {
	int r = 0;
	for(int i = 0; i < strlen(msg); i++)
	{
		int b1 = msg[i];		
		int b2 = extract_int8(c->allele,i*8);		
		for(size_t j = 0; j < 8; j++)
		{
			int bit1 = (b1 >> j) & 0x1;
			int bit2 = (b2 >> j) & 0x1;
			if(bit1 == bit2) r++;
		}
			
	}	
	return r;
}

double chromosome_to_double(struct chromosome *c)
{
	int i;
	double val;

	assert(c);

	for (i = 0, val = 0; i < chromosome_get_len(c); i++)
		if (chromosome_get_allele(c, i) == 1)
			val *= pow(2,i);

	return val;
}
void objective2(struct individual *i)
{
	double value=0;
	struct chromosome *c;
	Test * p;
	assert(i);	
	c = individual_get_chromosome(i);
	assert(c);	
	//value = extract_double(c->allele,0,64,-5,5);			
	p = (Test*)c->allele;
	value = p->x[0]*p->x[0] - p->x[1]*p->x[1];
	if(std::isinf(value)) value = 0;
	individual_set_fitness(i, value);
}

void test_ga2()
{
	Parameter<Test> p(100, 8*sizeof(Test), 500, 100, 0.995, 1e-6, 
			GA_S_ROULETTE_WHEEL, GA_X_SINGLE_POINT, objective2);
	p.set_report_strategy(GA_R_HUMAN_READABLE);
	p.evolve(10);
	//double james = extract_double(p.get_best_ever().individual().chromosome().p->allele,0,64,-5,5);
	printf("best %d=%f\n",p.p->best.generation,p.p->best.i->fitness);
	Test * x = p.decode_best();
	std::cout << x->x[0]*x->x[0] + x->x[1]*x->x[1] << std::endl;
	std::cout << p.get_best_ever().individual().chromosome().asString() << std::endl;
}

*/