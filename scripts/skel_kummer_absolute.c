read(strprintf("../src/functions_aux.c"));
read(strprintf("../src/functions_root_pol_gp.c"));
read(strprintf("../src/root_pol_multiemb.c"));
read(strprintf("../src/functions_exp.c"));

/* parameters to be defined */
/* exponent of kummer ext L, length of defining sequences, degree of eq. f(T),
   number of tests per field, version/type of eq : split or with only one root */

read(strprintf("../inputs/input_abs_kummer_%d_%d.c", EXPONENT, LENGTH));


VERSION = Str(VERSION);


default(parisizemax, 200000000000);
default(nbthreads, 1);

{
  ExperimentsKummer_LLL_abs(EXPONENT, LENGTH, DEGREE_EQ, vector_field_abs, NUMBER_TESTS, VERSION);
};

exit;
