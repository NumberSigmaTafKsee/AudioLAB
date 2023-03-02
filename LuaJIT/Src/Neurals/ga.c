#include "ga.h"


# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <limits.h>
#include <math.h>
#include <assert.h>
#include <sys/types.h>

#ifndef UINTBITS
# define UINTBITS	(sizeof(unsigned int) * CHAR_BIT)
#endif /* !UINTBITS */

#define NMEMB(len)		((int) ceil((double) len / UINTBITS))


extern void sgenrand(unsigned long seed);
extern double genrand(void);


void
random_seed(u_long seed)
{
	sgenrand(seed);
}


double
random_random(void)
{
	return genrand();
}


int
random_flip(float prob)
{
	return (random_random() <= prob) ? 1 : 0;
}


int
random_range(int low, int high)
{
	int i;

	if (low >= high)
		i = low;
	else {
		i = (random_random() * (high - low + 1)) + low;
		if (i > high)
			i = high;
	}

	return i;
}



/// chormosome

struct chromosome *
new_chromosome(size_t len)
{
	struct chromosome *c;

	if (!(c = malloc(sizeof(struct chromosome))))
		return NULL;

	if (!(c->allele = calloc(NMEMB(len), sizeof(int)))) {
		free(c);
		return NULL;
	}

	c->len = len;

	return c;
}

void
delete_chromosome(struct chromosome *self)
{
	assert(self);

	free(self->allele);
	free(self);
}


size_t
chromosome_get_len(struct chromosome *self)
{
	assert(self);

	return self->len;
}


void
set(unsigned int *i, unsigned int pos)
{
	*i |= (unsigned int) (1 << pos);
}

void
clear(unsigned int *i, size_t pos)
{
	*i &= (unsigned int) ~(1 << pos);
}

unsigned int
get(const unsigned int *i, size_t pos)
{
	return (unsigned int) ((*i >> pos) & 1); 
}

#define LOCUS(pos)		((int) floor((double) pos / UINTBITS))
#define BITOP(op, allele, pos)	\
	(op(allele[LOCUS(pos)], pos % UINTBITS))

void
chromosome_set_allele(struct chromosome *self, unsigned int pos)
{
	assert(self && pos < self->len);

	BITOP(set, &self->allele, pos);
}

void
chromosome_clear_allele(struct chromosome *self, unsigned int pos)
{
	assert(self && pos < self->len);

	BITOP(clear, &self->allele, pos);
}

unsigned int
chromosome_get_allele(struct chromosome *self, unsigned int pos)
{
	assert(self && pos < self->len);

	return BITOP(get, &self->allele, pos);
}

#undef LOCUS
#undef BITOP

void
chromosome_not_allele(struct chromosome *self, unsigned int pos)
{
	assert(self && pos < self->len);

	switch (chromosome_get_allele(self, pos)) {
	case 0:
		chromosome_set_allele(self, pos);
		break;
	case 1:
		chromosome_clear_allele(self, pos);
		break;
	}
}


char *
chromosome_as_string(struct chromosome *self)
{
	unsigned int i;
	char *s;

	assert(self);

	s = calloc(1, self->len + 1);
	if (!s)
		return NULL;

	for (i = 0; i < self->len; i++)
		s[i] = (chromosome_get_allele(self, i) == 1) ? '1' : '0';

	return s;
}


struct chromosome *
chromosome_dup(struct chromosome *self)
{
	struct chromosome *c;

	assert(self);

	c = new_chromosome(self->len);
	if (!c)
		return NULL;

	memcpy(c->allele, self->allele, sizeof(unsigned int) * NMEMB(self->len));

	return c;
}


void
chromosome_copy(struct chromosome *src, struct chromosome *dst,
		unsigned int src_pos, unsigned int dst_pos, size_t len)
{
	unsigned int i;
	size_t chrom_len;

	assert(src && dst);

	chrom_len = chromosome_get_len(src);
	assert(chrom_len == chromosome_get_len(dst));

	assert(src_pos < chrom_len);
	assert(dst_pos < chrom_len);
	assert(src_pos + len <= chrom_len);
	assert(dst_pos + len <= chrom_len);

	for (i = 0; i < len; i++) {
		switch (chromosome_get_allele(src, src_pos + i)) {
		case 0:
			chromosome_clear_allele(dst, dst_pos + i);
			break;
		case 1:
			chromosome_set_allele(dst, dst_pos + i);
			break;
		default:
			assert(0);	/* Just in case something went wrong */
			break;
		}
	}
}

/// individual

struct individual *
new_individual(size_t chrom_len)
{
	struct individual *i;

	i = calloc(1, sizeof(struct individual));
	if (!i)
		return NULL;

	i->chrom = new_chromosome(chrom_len);
	if (!i->chrom) {
		free(i);
		return NULL;
	}

	return i;
}

void
delete_individual(struct individual *self)
{
	assert(self);

	if (self->chrom) {
		/* See individual_dup for details on why this can happen. */
		delete_chromosome(self->chrom);
	}
	free(self);
}


struct chromosome *
individual_get_chromosome(struct individual *self)
{
	assert(self);

	return self->chrom;
}

void
individual_set_fitness(struct individual *self, double fitness)
{
	assert(self);

	self->fitness = fitness;
}

double
individual_get_fitness(struct individual *self)
{
	assert(self);

	return self->fitness;
}


void
individual_random(struct individual *self)
{
	unsigned int i;
	struct chromosome *c;
	void (*fn)(struct chromosome *, unsigned int);

	c = self->chrom;

	/* Slightly naive? Well... maybe */
	for (i = 0; i < chromosome_get_len(c); i++) {
		fn = (random_flip(0.5) == 1)
			? chromosome_set_allele
			: chromosome_clear_allele;

		fn(c, i);
	}
}


int
individual_compare(struct individual **lhs, struct individual **rhs)
{
	double f1, f2;
	int status = 0;

	f1 = individual_get_fitness(*lhs);
	f2 = individual_get_fitness(*rhs);
	
	if (f1 < f2)
		status = -1;
	else if (f1 == f2)
		status = 0;
	else if (f1 > f2)
		status = 1;

	return status;
}


struct individual *
individual_dup(struct individual *self)
{
	size_t chrom_len;
	struct individual *i;
	struct chromosome *c;

	assert(self);

	c = individual_get_chromosome(self);
	if (!c)
		return NULL;

	chrom_len = chromosome_get_len(c);

	i = new_individual(chrom_len);
	if (!i)
		return NULL;

	delete_chromosome(i->chrom);
	i->chrom = chromosome_dup(self->chrom);

	if (!i->chrom) {
		/*
		 * At this point we have an individual without a chromosome.
		 * That's the reason for the conditional in delete_individual.
		 */
		delete_individual(i);
		return NULL;
	}

	i->fitness = self->fitness;
		
	return i;
}


void
individual_print(struct individual *self, FILE *fp)
{
	char *s;

	assert(self && fp);

	s = chromosome_as_string(individual_get_chromosome(self));

	fprintf(fp, "%s (%.10f)\n", s, self->fitness);

	free(s);
}

/// population
struct population *
new_population(size_t len)
{
	struct population *p;

	assert(len > 0);

	p = calloc(1, sizeof(struct population));
	if (!p)
		return NULL;

	p->len = len;
	p->pop = calloc(p->len, sizeof(struct individual *));
	if (!p->pop) {
		free(p);
		return NULL;
	}

	return p;
}

void
delete_population(struct population *self)
{
	size_t i;

	assert(self);

	if (self->pop)
		for (i = 0; i < self->len; i++)
			if (self->pop[i]) delete_individual(self->pop[i]);
	free(self->pop);
	free(self->select_data);
	free(self);
}


struct fitness_stats *
population_get_fitness_stats(struct population *self)
{
	assert(self);

	return &self->stats;
}


struct individual *
population_get_fittest(struct population *self)
{
	assert(self);
	
	qsort(self->pop, self->len, sizeof(struct individual *),
		(int (*)(const void *, const void *)) individual_compare);
	return self->pop[self->len-1];	
}


void
population_compute_fitness_stats(struct population *self)
{
	int i;
	struct individual **pop;
	struct fitness_stats *stats;

	assert(self);

	stats = &self->stats;
	memset(stats, 0, sizeof(struct fitness_stats));

	pop = self->pop;

	/*
	 * This function assumes the population array has been sorted in
	 * advance (either by population_get_fittest or other means).
	 */
	stats->minimum = individual_get_fitness(pop[0]);
	stats->maximum = individual_get_fitness(pop[self->len - 1]);

	for (i = 0; i < self->len; i++)
		stats->total += individual_get_fitness(pop[i]);

	stats->average = stats->total / (double) self->len;
}


void
population_print(struct population *self, FILE *fp)
{
	int i;

        individual_print(population_get_fittest(self), stdout);
	/* for (i = 0; i < self->len; i++) */
	/* 	individual_print(self->pop[i], fp); */
}

/// operators

struct topbottom_pairing_data {
	unsigned int current;
};

void
preselect_topbottom_pairing(struct population *pop)
{
	struct topbottom_pairing_data *data;

	assert(pop);	/* new_pop can't be NULL */

	data = calloc(1, sizeof(struct topbottom_pairing_data));
	assert(data);	/* XXX - Need xmalloc */

	data->current = 0;

	assert(!pop->select_data);
	pop->select_data = data;
}

void
select_topbottom_pairing(struct population *pop,
			struct individual **dad, struct individual **mom)
{
	struct individual **pool;
	struct topbottom_pairing_data *data;

	assert(pop && mom && dad);

	data = (struct topbottom_pairing_data *) pop->select_data;

	pool = pop->pop;

	*dad = pool[data->current++];
	*mom = pool[data->current++];

	assert(*dad && *mom);
}



struct roulette_wheel_data {
	double total;
};

void
preselect_roulette_wheel(struct population *pop)
{
	struct roulette_wheel_data *data;

	assert(pop);

	data = calloc(1, sizeof(struct roulette_wheel_data));
	assert(data);

	data->total = population_get_fitness_stats(pop)->total;

	assert(!pop->select_data);
	pop->select_data = data;
}

void
select_roulette_wheel_parent(struct population *pop, struct individual **parent)
{
	unsigned int i;
	double sum, pick = random_random();
	struct individual *cur_indiv;
	struct roulette_wheel_data *data;

	assert(pop && parent);

	data = (struct roulette_wheel_data *) pop->select_data;
	assert(data);

	for (i = 0, sum = 0.0, cur_indiv = NULL; (sum < pick) && (i < pop->len); i++) {
		cur_indiv = pop->pop[i];

		sum += individual_get_fitness(cur_indiv) / data->total;
	}

	*parent = cur_indiv;
}

void
select_roulette_wheel(struct population *pop,
			struct individual **dad, struct individual **mom)
{
	assert(pop && dad && mom);

	select_roulette_wheel_parent(pop, dad);
	select_roulette_wheel_parent(pop, mom);
}



struct tournament_data {
	struct individual *mom,*dad;
};

void
preselect_tournament(struct population *pop)
{
	
	
}

void
select_tournament(struct population *pop,
			struct individual **dad, struct individual **mom)
{
	/*
	assert(pop && dad && mom);
	struct tournament_data candidate[10];

	for(size_t i = 0; i < 10; i++)
		select_roulette_wheel(pop,&candidate[i]->dad, &candidate[i]->mom);
	int x = (rand() % 10]->dad;
	*dad = candidate[x]->dad;
	*mom = candidate(x]->mom;
	*/
}



void
crossover_single_point(struct individual *dad, struct individual *mom,
			struct individual **son, struct individual **daughter)
{
	unsigned int xsite;		/* Crossover point */
	size_t chrom_len;
	struct chromosome *cdad, *cmom, *cson, *cdaughter;

	assert(dad && mom && son && daughter);

	cdad = individual_get_chromosome(dad);
	cmom = individual_get_chromosome(mom);
	assert(cdad && cmom);

	chrom_len = chromosome_get_len(cdad);
	assert(chrom_len == chromosome_get_len(cdad));

	*son = new_individual(chrom_len);
	*daughter = new_individual(chrom_len);
	assert(*son && *daughter);

	cson = individual_get_chromosome(*son);
	cdaughter = individual_get_chromosome(*daughter);
	assert(cson && cdaughter);

	chromosome_copy(cdad, cson, 0, 0, chrom_len);
	chromosome_copy(cmom, cdaughter, 0, 0, chrom_len);

	xsite = random_range(0, chrom_len - 1);
	(*son)->xsite = (*daughter)->xsite = xsite;

	chromosome_copy(cmom, cson, xsite, xsite, chrom_len - xsite);
	chromosome_copy(cdad, cdaughter, xsite, xsite, chrom_len - xsite);
}

void
crossover_n_point(struct individual *dad, struct individual *mom,
			struct individual **son, struct individual **daughter)
{
		
}

void
crossover_uniform(struct individual *dad, struct individual *mom,
			struct individual **son, struct individual **daughter)
{

}

size_t 
mutate(float pmutation, struct individual *indiv)
{
	unsigned int i;
	size_t mutations;
	struct chromosome *c;

	assert(indiv);

	c = individual_get_chromosome(indiv);
	assert(c);

	for (i = 0, mutations = 0; i < chromosome_get_len(c); i++) {
		if (random_flip(pmutation) == 0)
			continue;

		chromosome_not_allele(c, i);
		++mutations;
	}

	return mutations;
}

/// ga


void ga_set_best_ever(struct ga *self, unsigned int generation,
				struct individual *candidate);

void ga_set_selection_strategy(struct ga *self,
			enum ga_selection_strategies selection_strategy);

void ga_set_crossover_strategy(struct ga *self,
			enum ga_crossover_strategies crossover_strategy);

void ga_cross(struct ga *self, struct individual *dad,
			struct individual *mom, struct individual **son,
			struct individual **daughter);

void ga_preselect(struct ga *self);

void ga_mutate(struct ga *self, struct individual *indiv);

void ga_report(struct ga *self, FILE *fp);


struct ga *
new_ga(unsigned int max_gen,
	size_t chrom_len,
	size_t initial, size_t normal,
	float pcrossover, float pmutation,
	enum ga_selection_strategies selection_strategy,
	enum ga_crossover_strategies crossover_strategy,
	objective_fn obj_fn)
{
	struct ga *g;

	assert(chrom_len > 0
		&& initial > 0 && normal > 0
		&& pcrossover > 0 && pmutation > 0
		&& obj_fn);

	if (initial < normal) {
		fprintf(stderr, "%s: The initial population must be at least "
			"as big as the normal population\n", __FUNCTION__);
		return NULL;
	}
	if ((initial % 2 != 0) || (normal % 2 != 0)) {
		fprintf(stderr, "%s: Populations must have an even number "
			"of individuals\n", __FUNCTION__);
		return NULL;
	}

	g = calloc(1, sizeof(struct ga));
	if (!g)
		return NULL;

	g->current = 1;
	g->max_gen = max_gen;

	g->chrom_len = chrom_len;

	g->initial = initial;
	g->normal = normal;

	g->pcrossover = pcrossover;
	g->pmutation = pmutation;

	g->crossovers = g->mutations = 0;

	ga_set_selection_strategy(g, selection_strategy);
	ga_set_crossover_strategy(g, crossover_strategy);
	g->mutate = mutate;

	g->old_pop = g->cur_pop = NULL;

	g->best.i = NULL;
	g->best.generation = 0;

	g->obj_fn = obj_fn;

	return g;
}

void
delete_ga(struct ga *self)
{
	assert(self);

	if (self->best.i)
		delete_individual(self->best.i);
	if (self->old_pop)
		delete_population(self->old_pop);
	if (self->cur_pop)
		delete_population(self->cur_pop);

	free(self);
}

void ga_step(struct ga *self)
{
	if(self->current == 0) ga_first(self);
	else ga_next(self);
	self->current++;
}

void
ga_evolve(struct ga *self, unsigned int max_gen)
{
	assert(self);

	ga_first(self);
	while (self->current <= self->max_gen)
		ga_next(self);
}


void
ga_first(struct ga *self)
{
	unsigned int i;
	struct population *cur_pop;
	struct individual *cur_indiv;

	assert(self && !self->cur_pop);

	cur_pop = self->cur_pop = new_population(self->initial);

	/*
	 * Generate a random population.
	 */
	for (cur_indiv = NULL, i = 0; i < self->cur_pop->len; i++) {
		cur_indiv = cur_pop->pop[i] = new_individual(self->chrom_len);

		individual_random(cur_indiv);

		self->obj_fn(cur_indiv);
	}

	ga_set_best_ever(self, self->current,
			population_get_fittest(cur_pop));

	population_compute_fitness_stats(cur_pop);

	ga_report(self, stdout);

	++self->current;
}

void
swap_populations(struct population **old, struct population *new)
{
	if (*old)
		delete_population(*old);
	*old = new;
}

void
ga_next(struct ga *self)
{
	unsigned int i;
	struct population *cur_pop;
	struct individual *dad, *mom, *son, *daughter;

	assert(self && self->cur_pop);

	swap_populations(&self->old_pop, self->cur_pop);

	cur_pop = self->cur_pop = new_population(self->normal);

	ga_preselect(self);

	for (i = 0; i < self->normal; ) {
		self->select(self->old_pop, &dad, &mom);

		ga_cross(self, dad, mom, &son, &daughter);

		ga_mutate(self, son);
		ga_mutate(self, daughter);

		cur_pop->pop[i++] = son;
		cur_pop->pop[i++] = daughter;

		self->obj_fn(son);
		self->obj_fn(daughter);
	};

	ga_set_best_ever(self, self->current,
			population_get_fittest(cur_pop));

	population_compute_fitness_stats(cur_pop);

	ga_report(self, stdout);

	++self->current;
}


void
ga_preselect(struct ga *self)
{
	if (self->preselect)
		self->preselect(self->old_pop);
}

void
ga_cross(struct ga *self, struct individual *dad, struct individual *mom,
	struct individual **son, struct individual **daughter)
{
	if (random_flip(self->pcrossover) == 1) {
		self->cross(dad, mom, son, daughter);
		++self->crossovers;
	} else {
		*son = individual_dup(dad);
		*daughter = individual_dup(mom);
	}
}

void
ga_mutate(struct ga *self, struct individual *indiv)
{
	self->mutations += self->mutate(self->pmutation, indiv);
}


void
ga_set_best_ever(struct ga *self, unsigned int generation,
		struct individual *candidate)
{
	int status, should_change = 1;
	struct individual *best;

	assert(self && candidate);

	best = self->best.i;
	if (best) {
		if(candidate->fitness > best->fitness) status=1;
		else status=0;
	}
	
	/*
	 * In case the candidate is actually better than the best individual so
	 * far (or in case there's no best yet)...
	 */
	if (should_change) {
		if (best)			
			delete_individual(best);
		
		self->best.i = individual_dup(candidate);
		self->best.generation = generation;
	}
}

struct fittest *
ga_get_best_ever(struct ga *self)
{
	assert(self);

	return &self->best;
}


void
ga_set_selection_strategy(struct ga *self,
			enum ga_selection_strategies selection_strategy)
{
	preselection_fn preselect = NULL;
	selection_fn select = NULL;

	assert(self);

	switch (selection_strategy) {
	case GA_S_TOPBOTTOM_PAIRING:
		preselect = preselect_topbottom_pairing;
		select = select_topbottom_pairing;
		break;
	case GA_S_ROULETTE_WHEEL:
		preselect = preselect_roulette_wheel;
		select = select_roulette_wheel;
		break;
	case GA_S_TOURNAMENT:
		preselect = preselect_tournament;
		select = select_tournament;
		break;
	}

	self->preselect = preselect;
	self->select = select;
}

void
ga_set_crossover_strategy(struct ga *self,
			enum ga_crossover_strategies crossover_strategy)
{
	crossover_fn cross = NULL;

	assert(self);

	switch (crossover_strategy) {
	case GA_X_SINGLE_POINT:
		cross = crossover_single_point;
		break;
	case GA_X_N_POINT:
		cross = crossover_n_point;
		break;
	case GA_X_UNIFORM:
		cross = crossover_uniform;
		break;
	}

	self->cross = cross;
}


void
ga_report(struct ga *self, FILE *fp)
{
	assert(self && fp);

	if (self->report)
		self->report(self, fp);
}


void report_human_readable(struct ga *self, FILE *fp);
void report_graph(struct ga *self, FILE *fp);


void
ga_set_report_strategy(struct ga *self,
			enum ga_report_strategies report_strategy)
{
	report_fn report = NULL;

	assert(self);

	switch (report_strategy) {
	case GA_R_NOREPORT:
		break;
	case GA_R_HUMAN_READABLE:
		report = report_human_readable;
		break;
	case GA_R_GRAPH:
		report = report_graph;
		break;
	}

	self->report = report;
}


void
report_human_readable(struct ga *self, FILE *fp)
{
	struct fitness_stats *stats;

	stats = population_get_fitness_stats(self->cur_pop);
	assert(stats);

	fprintf(fp, ".----------------------------------------------------------------\n"); 
	fprintf(fp, "| Generation number: %3u / %3u\n", self->current, self->max_gen);
	fprintf(fp, "| Crossovers: %3u\tMutations: %3u\n", self->crossovers, self->mutations);
	fprintf(fp, "| Fittest individuals:\n");
	fprintf(fp, "|   so far (%3u):\t", self->best.generation);
	individual_print(self->best.i, fp);
	fprintf(fp, "|   in this generation:\t");
	individual_print(population_get_fittest(self->cur_pop), fp);
	fprintf(fp, "| Minimum fitness:\t%f\n", stats->minimum);
	fprintf(fp, "| Average fitness:\t%f\n", stats->average);
	fprintf(fp, "| Maximum fitness:\t%f\n", stats->maximum);
	fprintf(fp, "| Total fitness:\t%f\n", stats->total);
	fprintf(fp, "`----------------------------------------------------------------\n"); 
}

void
report_graph(struct ga *self, FILE *fp)
{
	struct fitness_stats *stats;

	stats = population_get_fitness_stats(self->cur_pop);
	assert(stats);

	fprintf(fp, "%u %f\n", self->current - 1, stats->average);
}

