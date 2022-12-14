read(strprintf("../src/functions_aux.c"));
read(strprintf("../src/functions_root_pol_gp.c"));
read(strprintf("../src/root_pol_multiemb.c"));
read(strprintf("../src/functions_exp.c"));



/* parameters to be defined */
/* dimensions dims, size pols, size of roots, number of tests */

/* size pols is chosen to be [1,1] */
SIZE_POL_VEC = [1, 1];
DIM_VEC = [DIM_G, DIM_E];

default(parisizemax, 200000000000);
default(nbthreads, 1);

{
  ExperimentsRelativeCompar_sizeroots(DIM_VEC, SIZE_POL_VEC, DEGREE_EQ, NUMBER_TESTS);
};

exit;
