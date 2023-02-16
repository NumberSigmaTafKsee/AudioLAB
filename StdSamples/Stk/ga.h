#pragma once

/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Name of package */
#define PACKAGE "libeve"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "jmbr@superadditive.com"

/* Define to the full name of this package. */
#define PACKAGE_NAME "libeve"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "libeve 0.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "libeve"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.1"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "0.1"

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

struct chromosome {
	size_t len;			/* Length (in bits) */
	unsigned int *allele;		/* Genetic material */
};


struct individual {
	struct chromosome *chrom;	/* Chromosome string */
	double fitness;			/* Fitness for this individual */
	unsigned int xsite;		/* Crossover site */
};



/*
 * Records the best individual so far and remembers the generation where it
 * appeared.
 */
struct fittest {
	struct individual *i;
	unsigned generation;
};


/*
 * Statistics about fitness in a generation.
 */
struct fitness_stats {
	double total;		/* Total fitness of a population */
	double minimum;		/* Less fit individual's fitness */
	double average;		/* Average fitness per individual */
	double maximum;		/* Best individual's fitness */
};

struct population {
	size_t len;			/* Population size */
	struct individual **pop;	/* Array of individuals */
	struct individual *fittest;	/* Fittest individual */
	struct fitness_stats stats;	/* Fitness statistics */
	size_t mutations;		/* Number of mutations */
	size_t crossovers;		/* Number of crossovers */

	void *select_data;		/* Selector-specific data */
};


typedef void (*preselection_fn)(struct population *);

typedef void (*selection_fn)(struct population *,
				struct individual **, struct individual **);

typedef void (*crossover_fn)(struct individual *, struct individual *,
				struct individual **, struct individual **);

typedef size_t (*mutation_fn)(float, struct individual *);

typedef void (*report_fn)(struct ga *, FILE *);

typedef void (*objective_fn)(struct individual *);

/*
 * A genetic algorithm.
 */
struct ga {
	unsigned int current;		/* Current generation index */
	unsigned int max_gen;		/* Maximum generation number */

	size_t chrom_len;		/* Length of a chromosome */

	size_t initial;			/* Individuals at 1st generation */
	size_t normal;			/* Individuals in a generation */

	float pcrossover;		/* Probability of crossover */
	float pmutation;		/* Probability of mutation */
	
	size_t crossovers;		/* Ammount of crossovers */
	size_t mutations;		/* Ammount of mutations */

	struct population *old_pop;	/* Previous population */
	struct population *cur_pop;	/* Current population */

	struct fittest best;		/* Best individual so far */

	objective_fn obj_fn;		/* User-defined objective function */

	preselection_fn preselect;	/* Preselection operator */
	selection_fn select;		/* Selection operator */

	crossover_fn cross;		/* Crossover operator */

	mutation_fn mutate;		/* Mutation operator */

	report_fn report;		/* Concrete report strategy */
};


enum ga_selection_strategies {
	GA_S_TOPBOTTOM_PAIRING,
	GA_S_ROULETTE_WHEEL,
	GA_S_TOURNAMENT
};

enum ga_crossover_strategies {
	GA_X_SINGLE_POINT,
	GA_X_N_POINT,
	GA_X_UNIFORM
};

enum ga_report_strategies {
	GA_R_NOREPORT = 0,
	GA_R_HUMAN_READABLE,
	GA_R_GRAPH
};



#ifdef __cplusplus
extern "C" {
#endif

struct ga *new_ga(unsigned int max_gen,
			size_t chrom_len,
			size_t initial, size_t normal,
			float pcrossover, float pmutation,
			enum ga_selection_strategies selection_strategy,
			enum ga_crossover_strategies crossover_strategy,
			objective_fn obj_fn);
void delete_ga(struct ga *self);

void ga_set_report_strategy(struct ga *self, enum ga_report_strategies report_strategy);

void ga_first(struct ga *self);
void ga_next(struct ga *self);
void ga_evolve(struct ga *self, unsigned maxgen);

struct fittest *ga_get_best_ever(struct ga *self);




void preselect_topbottom_pairing(struct population *pop);
void select_topbottom_pairing(struct population *pop,
					struct individual **dad,
					struct individual **mom);

void preselect_roulette_wheel(struct population *pop);
void select_roulette_wheel(struct population *pop,
					struct individual **dad,
					struct individual **mom);

void preselect_tournament(struct population *pop);
void select_tournament(struct population *pop,
					struct individual **dad,
					struct individual **mom);


void crossover_single_point(struct individual *dad,
					struct individual *mom,
					struct individual **son,
					struct individual **daughter);
void crossover_n_point(struct individual *dad, struct individual *mom,
				struct individual **son,
				struct individual **daughter);
void crossover_uniform(struct individual *dad, struct individual *mom,
				struct individual **son,
				struct individual **daughter);


size_t mutate(float pmutation, struct individual *indiv);



struct chromosome *chromosome_dup(struct chromosome *self);

void chromosome_copy(struct chromosome *src, struct chromosome *dst,
				unsigned int src_pos, unsigned int dst_pos, size_t len);



struct chromosome *new_chromosome(size_t len);
void delete_chromosome(struct chromosome *self);

size_t chromosome_get_len(struct chromosome *self);

void chromosome_set_allele(struct chromosome *self, unsigned int pos);
void chromosome_clear_allele(struct chromosome *self, unsigned int pos);
void chromosome_not_allele(struct chromosome *self, unsigned int pos);
unsigned int chromosome_get_allele(struct chromosome *self, unsigned int pos);

char *chromosome_as_string(struct chromosome *self);

/* Generates an individual with random genetic material */
void individual_random(struct individual *self);

int individual_compare(struct individual **lhs, struct individual **rhs);

struct individual *individual_dup(struct individual *self);

void individual_print(struct individual *self, FILE *fp);

struct individual *new_individual(size_t chrom_len);
void delete_individual(struct individual *self);

struct chromosome *individual_get_chromosome(struct individual *self);

void individual_set_fitness(struct individual *self, double fitness);
double individual_get_fitness(struct individual *self);




void population_compute_fitness_stats(struct population *self);

struct individual *population_get_fittest(struct population *self);

void population_print(struct population *self, FILE *fp);

void random_seed(u_long seed);

double random_random(void);

int random_flip(float prob);

int random_range(int low, int high);

#ifdef __cplusplus
}
#endif
