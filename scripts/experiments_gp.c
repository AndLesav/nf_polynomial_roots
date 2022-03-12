\r functions_aux.c;
\r functions_root_pol_gp.c;
\r root_pol_multiemb.c;
\r functions_exp.c;

/* \r input_rel_kummer.c */
/* \r input_global_kummer.c */


/* NUMBER_FIELDS = 10; */
NUMBER_TESTS = 20;

SIZE_POL_FIELD = 1;
DIM = 125;

SIZE_ROOT = 50;
FRAC_ROOTS = 25;

TYPE_FIELD = "complex";

TYPE_EQ = "";

SIZE_POL_VEC = [1, 1];


DIM_VEC = [64, 2];

/* DEGREE_POL = DIM_VEC[2]; */
DEGREE_POL = 50;

VERSION = "single";

/* /\* read files of precomputed polynomials defining number fields *\/ */
/* read(strprintf("input_rel_%d_%d.c", DIM_VEC[1], DIM_VEC[2])); */
/* read(strprintf("input_abs_%d_%d.c", DIM_VEC[1], DIM_VEC[2])); */

/* /\* read files of precomputed polynomials defining number fields *\/ */
/* read(strprintf("input_rel_reldeg_%d.c", DIM_VEC[1])); */
/* read(strprintf("input_abs_reldeg_%d.c", DIM_VEC[1])); */

/* read files of precomputed polynomials defining number fields -- kummer fields */
read(strprintf("input_rel_kummer_%d.c", DIM_VEC[2]));
read(strprintf("input_abs_kummer_%d.c", DIM_VEC[2]));

default(parisizemax, 200000000000);
default(nbthreads, 1);

{

  /* ##################################################################### */
  /* ####################  experiments -- global  ###################### */
  /* ##################################################################### */
  
  my(a);
  /* a = ExperimentsPolField(15, SIZE_POL_FIELD, SIZE_ROOT, DEGREE_POL, NUMBER_TESTS, TYPE_FIELD); */
  /* a = ExperimentsSizeRoots(DIM, SIZE_POL_FIELD, DEGREE_POL, NUMBER_TESTS, TYPE_FIELD); */
  /* a = ExperimentsEquation(DIM, SIZE_POL_FIELD, SIZE_ROOT, FRAC_ROOTS, NUMBER_TESTS, TYPE_FIELD); */
  /* a = ExperimentsNumberSolutions(DIM, SIZE_POL_FIELD, SIZE_ROOT, FRAC_ROOTS, NUMBER_TESTS, TYPE_FIELD, TYPE_EQ); */
  

  
  /* ####################################################################### */
  /* ######################  experiments -- relative  ######################### */
  /* ######################################################################### */

  /* COMPARISON OF RELATIVE METHODS  */

  /* a = ExperimentsRelativeCompar_deg(DIM_VEC, SIZE_POL_VEC, SIZE_ROOT, NUMBER_TESTS); */
  /* a = ExperimentsRelativeCompar_sizeroots(DIM_VEC, SIZE_POL_VEC, DEGREE_POL, NUMBER_TESTS); */

  /* **************************************************************************** */

  /* COMPARISON GLOBAL / RELATIVE METHODS  -- GENERIC CASE */

  /* equation */
  /* a = ExperimentsEquation_rel_abs(DIM_VEC, vector_field_rel, vector_field_abs, DEGREE_POL, SIZE_ROOT, NUMBER_TESTS); */


  /* relative degree */
  /* a = ExperimentsRelDegree_rel_abs(DIM_VEC[1], DEGREE_POL, vector_field_rel, vector_field_abs, SIZE_ROOT, NUMBER_TESTS); */

  
  /* **************************************************************************** */

  /* COMPARISON GLOBAL / RELATIVE METHODS +  NFROOTS -- KUMMER (NFROOTS is "bad")  */

  /* LLL abs version */
  /* ExperimentsKummer_LLL_abs(DIM_VEC, DEGREE_POL, vector_field_abs,  NUMBER_TESTS, VERSION); */

  /* LLL rel */
  /* ExperimentsKummer_LLL_rel(DIM_VEC, DEGREE_POL, vector_field_rel,  NUMBER_TESTS, VERSION); */

  /* gp version */
  /* ExperimentsKummer_GP(DIM_VEC, DEGREE_POL, vector_field_abs, NUMBER_TESTS, VERSION); */
  
  
  /* a = ExperimentsEquation_kummer(DIM_VEC, DEGREE_POL, vector_field_rel, SIZE_ROOT, NUMBER_TESTS); */
  /* a = ExperimentsEquation_kummer_abs(DIM_VEC, DEGREE_POL, vector_field_abs, SIZE_ROOT, NUMBER_TESTS); */



  /* ############################################################################## */
  /* ############################################################################## */


  /* OLD STUFF */

  /* a = ExperimentsRoots_kummer(DIM_VEC, DEGREE_POL, vector_field_rel, NUMBER_TESTS); */
  /* a = ExperimentsRoots_kummer_gp(DIM_VEC, DEGREE_POL, vector_field_abs, NUMBER_TESTS); */
  /* a = ExperimentsRoots_kummer_abs(DIM_VEC, DEGREE_POL, vector_field_abs, NUMBER_TESTS); */



  /* ******************************************************** */

  /* COMPARISON GLOBAL / RELATIVE METHODS + NFROOTS --  CYCLOTOMICS (NFROOTS is "good")  */
  ExperimentsCyclotomics(DIM, DEGREE_POL, SIZE_ROOT, NUMBER_TESTS);
  
  
};

exit;
