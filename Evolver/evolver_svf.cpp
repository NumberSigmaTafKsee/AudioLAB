#include <vector>
#include <complex>
#include "evolver.hpp"
#include "Plot.hpp"



#include "Analog/VAAnalogSVF.hpp"
#include "Analog/VABlitOscillators.hpp"

#define LOOP(start,end) for(size_t i = start; i < end; i++)

Default noise;
Plot<float> plot;

Analog::Filters::AnalogSVF::AnalogSVF filter(44100,12500,2.5);
Analog::Oscillators::Blit::BlitSaw osc;
std::vector<floatType> example,train;

struct SVF
{
	uint64_t cut;
	uint64_t res;
};

std::vector<floatType> genFilter(struct chromosome *c)
{
	SVF * p = (SVF*)c->allele;
	floatType cutoff = (floatType)p->cut/(floatType)UINT64_MAX;
	floatType q	     = (floatType)p->res/(floatType)UINT64_MAX;
	filter.reset();
    filter.setCutoff(cutoff*22050);
    filter.setQ(q*10.0);
    std::cout << cutoff*22050 << "," << q*10.0 << std::endl;
	std::vector<floatType> out(256);	
	LOOP(0,256) out[i] = filter.Tick(example[i]);
	return out;
}

int isFirst=0;
void objective(struct individual *i)
{
	double value=0;
	struct chromosome *c;	
	assert(i);	
	c = individual_get_chromosome(i);
	assert(c);		
	
	std::vector<floatType> out = genFilter(c);
	for(size_t i = 0; i < 256; i++)
	{					
		value += fabs(train[i] - out[i]);
	}
	
 	individual_set_fitness(i, -value);
 	
}

void test_svf()
{    
	plot.setstyle("lines");
	example.resize(256);
	train.resize(256);
	
    for(size_t i = 0; i < 256; i++) {
		example[i] = osc.Tick();
		train[i]   = filter.Tick(example[i]);
	}

		
    plot.plot_x(example.data(),example.size(),"examples");
	plot.plot_x(train.data(),train.size(),"training");
	
    GA ga(250, sizeof(SVF)*8, 500, 500, 0.9, 1e-6,GA_S_ROULETTE_WHEEL, GA_X_SINGLE_POINT, objective);
	ga.set_report_strategy(GA_R_HUMAN_READABLE);
	ga.evolve(10);	
	
	GAFittest fit = ga.best();
	GAIndividual ind = fit.individual();
	GAChromosome chrm= ind.chromosome();
	std::vector<floatType> out = genFilter(chrm.p);
	plot.plot_x(out.data(),out.size(),"genetic");
	sleep(1000);
}

int
main(int argc, char *argv[])
{
	srand(time(NULL));
	test_svf();
	
	exit(EXIT_SUCCESS);
}
