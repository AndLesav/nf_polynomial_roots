/* read(strprintf("../src/functions_aux.c")); */
read(strprintf("../src/functions_root_pol_gp.c"));
/* read(strprintf("../src/root_pol_multiemb.c")); */
read(strprintf("../src/functions_exp.c"));

/* parameters to be defined */
/* dimension start [K:QQ], size of defining pol P_K(X), size of roots, 
   degree of eq. f(T),   number of tests per field, type of field (r1 > 0) */

/* type of eq is assumed to be : split */

TYPE_FIELD = Str(TYPE_FIELD);
TYPE_EQ = Str(TYPE_EQ);

default(parisizemax, 200000000000);
default(nbthreads, 1);

{
  ExperimentsPolField(DIM_START, SIZE_POL_FIELD, SIZE_ROOTS, DEGREE_EQ, NUMBER_TESTS,\
		      NUMBER_DIM, TYPE_FIELD);
};

exit;
