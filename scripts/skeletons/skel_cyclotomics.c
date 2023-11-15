read(strprintf("../src/functions_aux.c"));
read(strprintf("../src/functions_root_pol_gp.c"));
/* read(strprintf("../src/root_pol_multiemb.c")); */
read(strprintf("../src/functions_exp.c"));

/* parameters to be defined :
   dimension max [K:QQ], size of roots, degree of eq. f(T), number of tests per field,
   interval of conductors.
   we assume split equations */

/* default(parisizemax, 200000000000); */
default(nbthreads, 1);

CONDS = [CONDS_MIN, CONDS_MAX];

{
  ExperimentsCyclotomicsNew(DIM_MAX, DEGREE_EQ, SIZE_ROOTS, NUMBER_TESTS, CONDS, "split");
};

exit;
