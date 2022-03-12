read(strprintf("../src/functions_aux.c"));
read(strprintf("../src/functions_root_pol_gp.c"));
read(strprintf("../src/root_pol_multiemb.c"));
read(strprintf("../src/functions_exp.c"));



/* parameters to be defined */
/* dimension [K:QQ], size of defining pol P_K(X), degree of eq. f(T),
   number of tests per field, type of field (r1 > 0) */

/* type of eq is assumed to be : split */

TYPE_FIELD = Str(TYPE_FIELD);
TYPE_EQ = Str(TYPE_EQ);


default(parisizemax, 200000000000);
default(nbthreads, 1);

{
  ExperimentsNumberSolutions(DIM, TYPE_FIELD, FRAC_ROOTS, SIZE_POL_FIELD, SIZE_ROOTS, \
			     NUMBER_TESTS, TYPE_EQ);
};

exit;
