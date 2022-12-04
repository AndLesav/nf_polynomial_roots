load "../src/functions_root_pol.m";
load "../src/general.m";
load "../src/functions_exp.m";


/* parameters to be defined */
/* rel dim of ext L/K, min and max dim of K, degree of eq. f(T),
   number of tests per field, version/type of eq : split or with only one root */

/* #################################################### */
/* #################################################### */
/* #################################################### */

SetClassGroupBounds("GRH");
p<z> := PolynomialRing(Integers());
q<y> := PolynomialRing(Integers());

if VERSION eq 0 then
    VERSION := "single";
elif VERSION eq 1 then
    VERSION := "split";
end if;

s := Sprintf("../inputs/input_reldeg_fixed_%o_%o.m", DIM_EXT, DIM_G_M);
print(s);
a := Read(s);
vector_field := eval a;


ExperimentsRelDegreeFixed_Magma(DIM_EXT, DIM_G_m, DIM_G_M, DEGREE_EQ,
				vector_field: number_tests:=NUMBER_TESTS,
					      max_size_roots := SIZE_ROOTS, 
					      version:=VERSION);

exit;
