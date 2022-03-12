p<z> := PolynomialRing(Integers());
q<y> := PolynomialRing(p);

print p, q;


NUMBER_FIELDS := 5*50;

load "functions_root_pol.m";


vector_field := [];
list_field_rel := [ ];

GROUND_DIMS := [15, 15, 15, 15];
EXT_DIMS := [2, 3, 4 ,5];


for ind in [1..#EXT_DIMS] do

    GROUND_DIM := GROUND_DIMS[ind];
    EXT_DIM := EXT_DIMS[ind];
    
    string_input_abs := Sprintf("input_abs_%o_%o.c", GROUND_DIM, EXT_DIM);;
    string_input_rel := Sprintf("input_rel_%o_%o.c", GROUND_DIM, EXT_DIM);;

    FILE := Open(string_input_abs, "w");
    fprintf FILE, "vector_field_abs = %o;\n", [];
    delete FILE;

    FILE := Open(string_input_rel, "w");
    fprintf FILE, "vector_field_rel = %o;\n", [];
    delete FILE;
        
    for j in [1..NUMBER_FIELDS] do

	print j,  " th number field over ", NUMBER_FIELDS;
	
	/* representation is a bit sparse for efficiency reasons */
	L, K := Pol_field_creation_tower([GROUND_DIM, EXT_DIM], 1 : supp := [1/5,1/2,1/2]);
	
	K := AbsoluteField(K);
	P := p!DefiningPolynomial(K);

	p<z> := PolynomialRing(Integers());
	q<y> := PolynomialRing(p);
	
	A := [ Eltseq(L[i]) : i in [1..#L] ];
	a := p!A[1];
	B := A[2];
	B := [p!Eltseq(B[i]) : i in [1..#B]];
	b := q!B;
	
	print a, b;
	
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
