/* read(strprintf("../src/functions_aux.c")); */
read(strprintf("../src/functions_root_pol_gp.c"));
/* read(strprintf("../src/root_pol_multiemb.c")); */
read(strprintf("../src/functions_exp.c"));

/* parameters to be defined */
/* dimension start [K:QQ], size of defining pol P_K(X), number of tests, 
   degree of eq. f(T), type of field (r1 > 0), number of dim */

/* type of eq is assumed to be : split */

TYPE_FIELD = Str(TYPE_FIELD);

default(parisizemax, 200000000000);
default(nbthreads, 1);

{
  ExperimentsProportionBadFields(DIM_START, SIZE_POL_FIELD, NUMBER_TESTS, \
				 DEGREE_EQ, TYPE_FIELD, NUMBER_DIM);
};

exit;
