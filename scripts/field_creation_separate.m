p<z> := PolynomialRing(Integers());
q<y> := PolynomialRing(p);

print p, q;

load "../src/functions_root_pol.m";
load "../src/general.m";


GROUND_DIM := DIM_G;
EXT_DIM := DIM_EXT;

vector_field_good := [];
vector_field_bad := [];

list_field_rel_good := [ ];
list_field_rel_bad := [ ];

string_input_abs_good := Sprintf("../inputs/input_separate_abs_good_%o_%o.c", GROUND_DIM, EXT_DIM);
string_input_rel_good := Sprintf("../inputs/input_separate_rel_good_%o_%o.c", GROUND_DIM, EXT_DIM);
string_input_abs_bad := Sprintf("../inputs/input_separate_abs_bad_%o_%o.c", GROUND_DIM, EXT_DIM);
string_input_rel_bad := Sprintf("../inputs/input_separate_rel_bad_%o_%o.c", GROUND_DIM, EXT_DIM);


FILE := Open(string_input_abs_good, "w");
fprintf FILE, "vector_field_abs = %o;\n", [];
delete FILE;
FILE := Open(string_input_abs_bad, "w");
fprintf FILE, "vector_field_abs = %o;\n", [];
delete FILE;

FILE := Open(string_input_rel_good, "w");
fprintf FILE, "vector_field_rel = %o;\n", [];
delete FILE;
FILE := Open(string_input_rel_bad, "w");
fprintf FILE, "vector_field_rel = %o;\n", [];
delete FILE;


for j in [1..NUMBER_FIELDS] do

    print j,  " th number field over ", NUMBER_FIELDS;

    /* representation is a bit sparse for efficiency reasons */
    L, K := Pol_field_creation_tower([GROUND_DIM, EXT_DIM], 5 : supp := [1/5,1/2,1/2]);

    K := AbsoluteField(K);
    P := p!DefiningPolynomial(K);

    bool, best_f := Is_good_field(K);
    
    p<z> := PolynomialRing(Integers());
    q<y> := PolynomialRing(p);
    
    A := [ Eltseq(L[i]) : i in [1..#L] ];
    a := p!A[1];
    B := A[2];
    B := [p!Eltseq(B[i]) : i in [1..#B]];
    b := q!B;
    
    print a, b;

    if bool then 
	Append(~vector_field_good, P);
	Append(~list_field_rel_good, [a, b]);
    else
	Append(~vector_field_bad, P);
	Append(~list_field_rel_bad, [a, b]);
    end if;


    if bool then
	string_input_abs := string_input_abs_good;
	string_input_rel := string_input_rel_good;
    else
	string_input_abs := string_input_abs_bad;
	string_input_rel := string_input_rel_bad;
    end if;

    FILE := Open(string_input_abs, "a");
    fprintf FILE, "vector_field_abs = concat(vector_field_abs, [%o, %o]); \n", q!P, best_f;
    delete FILE;
    
    FILE := Open(string_input_rel, "a");
    fprintf FILE, "vector_field_rel = concat(vector_field_rel, [[%o, %o], %o]); \n", p!a, q!b, best_f;
    delete FILE;
end for;
