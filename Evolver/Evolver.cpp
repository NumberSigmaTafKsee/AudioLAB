#include <cstddef>
#include <cmath>
#include <cstring>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <memory>

#define LOOP(i,start,end) for(size_t i = start; i < end; i++)

struct Random
{
	Random() {
		srand(time(NULL));
	}
		
	void seed(uint64_t seed)
	{
		srand(seed);
	}


	double random(void)
	{
		return (double)rand()/(double)RAND_MAX;
	}

	int randint(int lo, int up) {
		return (double)(up-lo)*random() + lo;
	}

	int flip(float prob)
	{
		return (random() <= prob) ? 1 : 0;
	}


	int range(int low, int high)
	{
		int i;

		if (low >= high)
			i = low;
		else {
			i = (random() * (high - low + 1)) + low;
			if (i > high)
				i = high;
		}

		return i;
	}
};

Random noise;


struct BitString
{
	std::vector<uint8_t> bits;
	int num_bits;
	
	BitString(int n) {
		num_bits = n;
		bits.resize(ceil((double)n/8));
		memset(bits.data(),0,bits.size());
	}
	size_t size() const { return num_bits; }
	void* data() { return (void*)bits.data(); }
	
	bool operator[](size_t pos) {
		int p   = pos/8;
		int x   = pos % 8;
		return (bits[p] >> x) & 0x1;
	}
	void clear(size_t pos) {
		int p   = pos/8;
		int x   = pos % 8;
		bits[p] &= ~(1 << x);
	}
	void set(size_t pos) {
		int p   = pos/8;
		int x   = pos % 8;		
		bits[p] |= 1 << x;
	}
	void flip(size_t pos) {
		int p   = pos/8;
		int x   = pos % 8;		
		if(bits[p] & (1 << x)) clear(pos);
		else set(pos);
	}
	void print() {
		for(size_t i = 0; i < num_bits; i++)
		{
			std::cout << (int)(*this)[i];
		}
		std::cout << std::endl;
	}
	void randomize() {
		for(size_t i = 0; i < num_bits; i++)
		{
			if(noise.flip(0.5)) set(i);
			else clear(i);
		}
	}
};

struct Object : public BitString
{
	double    fitness;
	uint64_t generation;
	
	Object(size_t n) : BitString(n) {
		fitness = 0;
		generation = 0;
	}	
	Object& operator = (const Object& o) {
		bits = o.bits;
		num_bits = o.num_bits;
		fitness = o.fitness;
		generation = o.generation;
		return *this;
	}
};
	
	
struct Population
{
	std::vector<std::shared_ptr<Object>> pop;
	double total;		/* Total fitness of a population */
	double minimum;		/* Less fit individual's fitness */
	double average;		/* Average fitness per individual */
	double maximum;		/* Best individual's fitness */
	size_t crossovers;
	size_t mutations;
	size_t len;
	
	Population(size_t n, size_t clen) {
		init(n,clen);
		total = 0;
		minimum = 0;
		maximum = 0;
		average = 0;
		crossovers = 0;
		mutations = 0;		
		len = clen;
	}
	void init(size_t n,size_t clen) {		
		pop.resize(n);
		for(size_t i = 0; i < n; i++)
			pop[i] = std::make_shared<Object>(clen);
	}
	Object* operator[](size_t i) {
		return pop[i].get();
	}
		
	size_t size() const { return pop.size(); }
	
	Object* select_parent();
	void 	select(Object** dad, Object** mom);
	
};
		
void crossover(Object & p1, Object & p2, Object & c1, Object & c2) 
{	
	size_t pos = noise.randint(1,p1.size());	
	
	for(size_t i = 0; i < pos; i++) 
	{		
		if(p1[i]) c2.set(i);
		else c2.clear(i);
		if(p2[i]) c1.set(i);
		else c1.clear(i);
	}	
}
void mutate(Object & c)
{
	size_t pos = noise.randint(0,c.size());	
	if(c[pos]) c.clear(pos);
	else c.set(pos);
}

Object* Population::select_parent()
{
	unsigned int i;
	double sum, pick = noise.random();
	Object* cur_indiv;
	
	for (i = 0, sum = 0.0, cur_indiv = NULL; (sum < pick) && (i < len); i++) {
		cur_indiv = (*this)[i];
		sum += cur_indiv->fitness / total;
	}
	return cur_indiv;
}

void Population::select(Object** dad, Object** mom)
{	
	*dad = select_parent();
	*mom = select_parent();
}

void swap(Population &a, Population & b)
{
	Population temp = a;
	a = b;
	b = temp;
}

struct ga {
	unsigned int current;		/* Current generation index */
	unsigned int max_gen;		/* Maximum generation number */
	bool show_report;
	
	size_t chrom_len;		/* Length of a chromosome */

	size_t initial;			/* Individuals at 1st generation */
	size_t normal;			/* Individuals in a generation */

	double pcrossover;		/* Probability of crossover */
	double pmutation;		/* Probability of mutation */
	
	size_t crossovers;		/* Ammount of crossovers */
	size_t mutations;		/* Ammount of mutations */

	std::shared_ptr<Population> old_pop;	/* Previous population */
	std::shared_ptr<Population> cur_pop;	/* Current population */

	std::shared_ptr<Object> best;		/* Best individual so far */

	std::function<void (Object&)> obj_fn;	
	
	ga(size_t max, size_t clen,size_t init, size_t norm, double cross, double mutations, std::function<double (Object&)> obj) {
		max_gen = max;
		current = 0;
		max_gen = 0;
		show_report = false;
		chrom_len = clen;
		initial = init;
		normal  = norm;
		obj_fn = obj;
		pcrossover = cross;
		pmutation = mutations;
		old_pop = std::make_shared<Population>(clen);
		cur_pop = std::make_shared<Population>(clen);
		best = std::make_shared<Object>(clen);
	}
	
	void Evolve() {
		First();
		for(current=1; current <= max_gen; current++)	
			Next();		
	}
	void First() {				
		Population *p;
		cur_pop.reset();
		cur_pop = std::make_shared<Population>(initial,chrom_len);
		p = cur_pop.get();	
		
	
		#pragma omp parallel for
		for (size_t i = 0; i < p->size(); i++) {
			Object* cur_indiv = (*p)[i];
			cur_indiv->randomize();		
			obj_fn(*cur_indiv);
		}
			
		if(cur_indiv->fitness > best->fitness)
			*best = *cur_indiv;
					
		p->compute_fitness_stats();
		report();

		++current;
	}
	void Next()
	{
		uint64_t i;
		Population *p,*old;
		Object *dad, *mom, *son, *daughter;
		
		swap(old_pop, cur_pop);
		cur_pop.reset();
		cur_pop = std::make_shared<Object>(normal,clen);
		old = old_pop.get();
		
		#pragma omp parallel for
		for (i = 0; i < normal; i+=2) {
			old_pop->select(dad, mom);
			son = std::make_shared<Object>(clen);
			daugher = = std::make_shared<Object>(clen);
			*son = *dad;
			*daughter = *mom;
			
			if(noise.random() < pcrossover)
				crossover(*dad, *mom, *son, *daughter);
	
			if(noise.random() < pmutation)
				mutation(*son);
				
			if(noise.random() < pmutation)
				mutation(*daughter);

			(*cur_pop)[i] 	= *son;
			(*cur_pop)[i+1] = *daughter;

			obj_fn(son);
			obj_fn(daughter);
		};

		*best = cur_pop->get_fittest();		
		cur_pop->compute_fitness_stats();
		report();

		++current;
	}

};	



int main()
{
	
}
