read(strprintf("../src/functions_root_pol_gp.c"));
read(strprintf("../src/functions_exp.c"));

/* parameters to be defined */
/* dimensions dims, size pols, size of roots, number of tests */

/* size pols is chosen to be [1,1] */
DIM_VEC = [DIM_G, DIM_E];

/* read files of precomputed polynomials defining number fields */
read(strprintf("../inputs/input_rel_nbsol_%d_%d.c", DIM_VEC[1], DIM_VEC[2]));
read(strprintf("../inputs/input_abs_nbsol_%d_%d.c", DIM_VEC[1], DIM_VEC[2]));

default(parisizemax, 200000000000);
default(nbthreads, 1);


/* CAREFUL:  NUMBER_TESTS NEEDS TO BE THE SAME AS THE ONE USED TO CREATE FIELDS */

{
  ExperimentsNbSolutions_abs_rel(DIM_VEC, vector_field_rel, vector_field_abs, DEGREE_EQ, SIZE_ROOTS, NUMBER_TESTS);
};

exit;
