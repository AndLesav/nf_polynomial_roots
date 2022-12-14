read(strprintf("../src/functions_root_pol_gp.c"));
read(strprintf("../src/functions_exp.c"));

/* parameters to be defined */
/* dimension max [K:QQ], size of roots, 
   degree of eq. f(T),   number of tests per field, prime, exp of extension */

default(parisizemax, 200000000000);
default(nbthreads, 1);


{
  ExperimentsPrimePowerCyclotomics(PRIME, DEGREE_EQ, SIZE_ROOTS, NUMBER_TESTS, \
				   DIM_MAX);
};

exit;
