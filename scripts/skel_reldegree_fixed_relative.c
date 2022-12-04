read(strprintf("../src/functions_aux.c"));
read(strprintf("../src/functions_root_pol_gp.c"));
read(strprintf("../src/root_pol_multiemb.c"));
read(strprintf("../src/functions_exp.c"));

/* parameters to be defined */
/* rel dim of ext L/K, min and max dim of K, degree of eq. f(T),
   number of tests per field, version/type of eq : split or with only one root */

read(strprintf("../inputs/input_rel_reldeg_fixed_%d_%d.c", DIM_EXT, DIM_G_M));

VERSION = Str(VERSION);


default(parisizemax, 200000000000);
default(nbthreads, 1);
{
  ExperimentsRelDegreeFixed_LLL_rel(DIM_EXT, DIM_G_m, DIM_G_M, DEGREE_EQ, \
				    vector_field_rel, NUMBER_TESTS, SIZE_ROOTS,\
				    VERSION);
};

exit;
