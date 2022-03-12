p<z> := PolynomialRing(Integers());
q<y> := PolynomialRing(p);

print p, q;

load "../src/functions_root_pol.m";

vector_field := [];
list_field_rel := [ ];

GROUND_DIMS := [DIM_G_m];
while #GROUND_DIMS lt NUMBER_DIM do
    Append(~GROUND_DIMS, GROUND_DIMS[#GROUND_DIMS]+10);
end while;

EXT_DIMS := [2, 3, 4, 5];

for GROUND_DIM in GROUND_DIMS do

    string_input_abs := Sprintf("../inputs/input_abs_reldeg_%o.c", GROUND_DIM);
    string_input_rel := Sprintf("../inputs/input_rel_reldeg_%o.c", GROUND_DIM);

    FILE := Open(string_input_abs, "w");
    fprintf FILE, "vector_field_abs = %o;\n", [];
    delete FILE;

    FILE := Open(string_input_rel, "w");
    fprintf FILE, "vector_field_rel = %o;\n", [];
    delete FILE;

    
    for EXT_DIM in EXT_DIMS do
	
	for j in [1..NUMBER_FIELDS] do
	    
	    print j,  " th number field over ", NUMBER_FIELDS;

	    /* representation is a bit sparse for efficiency reasons */
	    L, K := Pol_field_creation_tower([GROUND_DIM, EXT_DIM], 1 : supp := [1/3, 1/2, 1/2]);
	    
	    K := AbsoluteField(K);
	    P := p!DefiningPolynomial(K);

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
	    fprintf FILE, "vector_field_abs = concat(vector_field_abs, [%o]); \n", q!P;
	    delete FILE;
	    
	    FILE := Open(string_input_rel, "a");
	    fprintf FILE, "vector_field_rel = concat(vector_field_rel, [[%o, %o]]); \n", p!a, q!b;
	    delete FILE;

	end for;
    end for;
end for;
