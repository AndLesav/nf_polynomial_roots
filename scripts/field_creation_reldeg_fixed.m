p<z> := PolynomialRing(Integers());
q<y> := PolynomialRing(p);

print p, q;

load "../src/functions_root_pol.m";
load "../src/general.m";

vector_field := [];
list_field_rel := [ ];

EXT_DIM := DIM_EXT;
string_input_abs := Sprintf("../inputs/input_abs_reldeg_fixed_%o_%o.c", EXT_DIM, DIM_G_M);
string_input_rel := Sprintf("../inputs/input_rel_reldeg_fixed_%o_%o.c", EXT_DIM, DIM_G_M);
string_input_magma := Sprintf("../inputs/input_reldeg_fixed_%o_%o.m", EXT_DIM, DIM_G_M);

FILE := Open(string_input_abs, "w");
fprintf FILE, "vector_field_abs = %o;\n", [];
delete FILE;

FILE := Open(string_input_rel, "w");
fprintf FILE, "vector_field_rel = %o;\n", [];
delete FILE;

for i in [1..NUMBER_FIELDS] do

    print "\t", i,  " th number field over ", NUMBER_FIELDS;

    
    GROUND_DIM := Random(DIM_G_m, DIM_G_M);
    dim_abs := EXT_DIM * GROUND_DIM;
    
    /* representation can be chosen sparse for efficiency reasons */
    L, K := Pol_field_creation_tower([GROUND_DIM, EXT_DIM], 1 : supp := [1, 1, 1]);
    
    K := AbsoluteField(K);
    P := p!DefiningPolynomial(K);
    
    bool, best_f := Is_good_field(K: deg_eq := DEGREE_EQ);
    print(bool);
    if bool then
	bool := 1;
    else
	bool := 0;
    end if;
    
    p<z> := PolynomialRing(Integers());
    q<y> := PolynomialRing(p);
    
    A := [ Eltseq(L[i]) : i in [1..#L] ];
    a := p!A[1];
    B := A[2];
    B := [p!Eltseq(B[i]) : i in [1..#B]];
    b := q!B;
    
    Append(~vector_field, P);
    Append(~list_field_rel, [a, b]);
    
    FILE := Open(string_input_abs, "a");
    fprintf FILE, "vector_field_abs = concat(vector_field_abs, [%o, %o]); \n", q!P, bool;
    delete FILE;
    
    FILE := Open(string_input_rel, "a");
    fprintf FILE, "vector_field_rel = concat(vector_field_rel, [[%o, %o, %o]]); \n", p!a, q!b, bool;
    delete FILE;

end for;

FILE := Open(string_input_magma, "a");
fprintf FILE, "vector_field := %o; \n", vector_field;
fprintf FILE, "return vector_field;", vector_field;
delete FILE;
