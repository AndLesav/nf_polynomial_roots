p<z> := PolynomialRing(Integers());
q<y> := PolynomialRing(p);

print p, q;



load "../src/functions_root_pol.m";

vector_field := [];
list_field_rel := [ ];

exp := EXPONENT;
r := LENGTH;

GROUND_DIM := exp^(r-1);
EXT_DIM := exp;

string_input_abs := Sprintf("../inputs/input_abs_kummer_%o_%o.c", exp, r);
string_input_rel := Sprintf("../inputs/input_rel_kummer_%o_%o.c", exp, r);

FILE := Open(string_input_abs, "w");
fprintf FILE, "vector_field_abs = %o;\n", [];
delete FILE;

FILE := Open(string_input_rel, "w");
fprintf FILE, "vector_field_rel = %o;\n", [];
delete FILE;



for j in [1..NUMBER_FIELDS] do /* #SIZES * NUMBER_FIELDS */
    L, K := Pol_field_creation_kummer(exp, r, PRIMES_RANGE);
    
    K := AbsoluteField(K);
    P := p!DefiningPolynomial(K);
    p<z> := PolynomialRing(Integers());
    q<y> := PolynomialRing(p);
    

    A := [ L[i] : i in [1..#L-1] ];
    A := AbsoluteField(NumberField(A));
    a := p!DefiningPolynomial(A);
    
    b := q!L[#L];
        
    FILE := Open(string_input_abs, "a");
    fprintf FILE, "vector_field_abs = concat(vector_field_abs, [%o]); \n", q!P;
    delete FILE;
    
    FILE := Open(string_input_rel, "a");
    fprintf FILE, "vector_field_rel = concat(vector_field_rel, [[%o, subst(%o, z, y)]]); \n", p!a, q!b;
    delete FILE;

end for;

