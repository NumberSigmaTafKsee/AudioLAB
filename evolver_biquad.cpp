#include <vector>
#include <complex>
#include "evolver.hpp"
#include "Plot.hpp"
#include "ml_helpers.hpp"

Default noise;
Plot<float> plot;

Biquad biquad;
int FR_SIZE = 128;
std::vector<float> bq_example(FR_SIZE),bq_train(FR_SIZE);

struct Filter
{
    int64_t z[3];
    int64_t p[3];
    int64_t k;
};


std::vector<floatType> frequencyResponse() {
    const int ir_size = FR_SIZE;
    AudioDSP::FFTPlanRealFloat fftPlan(ir_size);
    std::vector<std::complex<floatType>> fr(ir_size);
    std::vector<floatType> r(ir_size);        
    biquad.clear();
    r[0] = biquad.Tick(1.0);
    for(size_t i = 1; i < ir_size; i++)
        r[i] = biquad.Tick(0);        
    AudioDSP::fft(fftPlan,r.data(),fr.data());  
    for(size_t i = 0; i < ir_size; i++)
    {
		fr[i] /= std::complex<floatType>(ir_size,ir_size);
	}  
    for(size_t i = 0; i < ir_size; i++) r[i] = std::abs(fr[i]);        
    return r;
}

void load(Filter * p) {
	floatType G = ((float)p->k / (float)INT64_MAX);
    for(size_t i = 0; i < 3; i++)	{
        biquad.z[i] = G*((float)p->z[i] / (float)INT64_MAX);
        biquad.p[i] = ((float)p->p[i] / (float)INT64_MAX);
        if(std::isnan(p->z[i]) || std::isinf(p->z[i])) p->z[i] = 0;
        if(std::isnan(p->p[i]) || std::isinf(p->z[i])) p->p[i] = 0;
    }	
}   
int isFirst=0;
void objective(struct individual *i)
{
	double value=0;
	struct chromosome *c;	
	assert(i);	
	c = individual_get_chromosome(i);
	assert(c);		
	Filter *p = (Filter*)c->allele;
	
	load(p);
    std::vector<floatType> fr = frequencyResponse();
    
    
	if(isFirst < 100) {
		isFirst++;
		value   = 0;
	}	
	else {
		for(size_t i = 0; i < FR_SIZE; i++)
		{					
			value += fabs(fr[i] - bq_train[i]) < 1e-6;
		}
	}
 	individual_set_fitness(i, value);
 	
}

void test_biquad()
{    
	plot.setstyle("lines");
    ButterworthHighpassFilter butter(2,44100.0);
    bq_example = bq_train = butter.frequency_response(FR_SIZE);
    
    plot.plot_x(bq_example.data(),bq_example.size(),"examples");
	plot.plot_x(bq_train.data(),bq_train.size(),"training");
	
    
    GA ga(250, sizeof(Filter)*8, 500, 500, 0.9, 1e-6,GA_S_ROULETTE_WHEEL, GA_X_SINGLE_POINT, objective);
	ga.set_report_strategy(GA_R_HUMAN_READABLE);
	ga.evolve(10);	
	Filter * p = (Filter*)(ga.p->best.i->chrom->allele);
	load(p);
	
	std::vector<floatType> out = frequencyResponse();	
	plot.plot_x(out.data(),out.size(),"output");
	
	sleep(1000);
}

int
main(int argc, char *argv[])
{
	srand(time(NULL));
	test_biquad();
	
	exit(EXIT_SUCCESS);
}
