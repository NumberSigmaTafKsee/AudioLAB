#include "evolver.hpp"
typedef float DspFloatType;
#include "Plot.hpp"
#include "stdsamples.hpp"
#include "ml_helpers.hpp"
#include "Luajit.hpp"

Default noise;
Plot<float> plot;
bool doSort = true;


///////////////////////////////////////////////////////////////
// Test area
///////////////////////////////////////////////////////////////
static void test_ga(void);
static void objective(struct individual *i);
static double chromosome_to_double(struct chromosome *c);
int check_bits(struct chromosome * c, const char * msg);


////////////////////
// Chebyshev
////////////////////

#define CHEBYS 8
struct Cheby
{
	int64_t A[CHEBYS];
};

int cmpf(const void * a, const void * b) {
	int64_t x = *(int64_t*)a;
	int64_t y = *(int64_t*)b;
	if(x < y) return 1;
	if(x > y) return -1;
	return 0;
}	

floatType chebyshev(floatType x, int64_t A[], int order)
{
	
	// To = 1
	// T1 = x
	// Tn = 2.x.Tn-1 - Tn-2
	// out = sum(Ai*Ti(x)) , i C {1,..,order}
	floatType Tn_2 = 1.0f;	
	floatType Tn_1 = x;
	floatType Tn;
	auto X  =   (((floatType)(A[0])/(floatType)INT64_MAX));
	floatType out = X*x;

	#pragma omp simd
	for(int n=2;n<=order;n++)
	{
		auto Q  =   (((floatType)(A[n-1])/(floatType)INT64_MAX));
		Tn    	=   2.0*x*Tn_1 - Tn_2;
		out    +=   Q*Tn;
		Tn_2 	=   Tn_1;
		Tn_1 	=   Tn;
	}
	
	return out;
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
	if(doSort) qsort(p->A,CHEBYS,sizeof(INT64_MAX),cmpf);
	#pragma omp simd	
 	for(size_t i = 0; i < 128; i++) {
		auto bobo = chebyshev(examples[i],p->A,CHEBYS) - chebyshev(0,p->A,CHEBYS);
		auto x = bobo; 
		out[i] = x;		
	}	
	
	RemoveDCNorm(out);
	
		
	#pragma omp simd
	for(size_t i = 0; i < 128; i++) {		
		auto x = fabs(out[i] - training[i]) < 0.1;
		value += x;
	}		
 	individual_set_fitness(i, value);
}

void test_cheby()
{	
	#pragma omp simd
    for(int i = 0; i < 128; i++) {
		auto sample = std::sin(2*M_PI*(double)i/128.0);
        examples.push_back(sample);
        training.push_back(FX::Distortion:: ArayaAndSuyama(sample));
    }
	plot.setstyle("lines");
    plot.plot_x(examples.data(),examples.size(),"examples");
	plot.plot_x(training.data(),training.size(),"training");
	GA ga(100, sizeof(Cheby)*8, 500, 50, 0.9, 1e-3,GA_S_ROULETTE_WHEEL, GA_X_SINGLE_POINT, objective);
	ga.set_report_strategy(GA_R_HUMAN_READABLE);
	ga.evolve(1000);
	Cheby * p = (Cheby*)(ga.p->best.i->chrom->allele);
	std::vector<float> out(128),temp;
	for(size_t i = 0; i < CHEBYS; i++) std::cout << (floatType)p->A[i]/(floatType)INT64_MAX << ",";
	std::cout << std::endl;
	#pragma omp simd
	for(size_t i = 0; i < 128; i++)
	{
		auto bobo = chebyshev(examples[i],p->A,CHEBYS) - chebyshev(0,p->A,CHEBYS);
		out[i] = bobo; //(bobo-1+0.5)/4;
	}
	
	RemoveDCNorm(out);
	
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
