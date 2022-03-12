/* read(strprintf("../src/functions_aux.c")); */
read(strprintf("../src/functions_root_pol_gp.c"));
/* read(strprintf("../src/root_pol_multiemb.c")); */
read(strprintf("../src/functions_exp.c"));

/* parameters to be defined */
/* dimensions dim of ground field, degree of equation, size of roots, number of tests */

/* read files of precomputed polynomials defining number fields */
read(strprintf("../inputs/input_rel_reldeg_%d.c", DIM_G));
read(strprintf("../inputs/input_abs_reldeg_%d.c", DIM_G));

default(parisizemax, 200000000000);
default(nbthreads, 1);

/* CAREFUL:  NUMBER_TESTS NEEDS TO BE THE SAME AS THE ONE USED TO CREATE FIELDS */


{
  ExperimentsRelDegree_abs_rel(DIM_G, DEGREE_EQ, vector_field_rel, vector_field_abs, SIZE_ROOTS, NUMBER_TESTS);
};

exit;
