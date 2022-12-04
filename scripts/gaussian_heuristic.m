load "../src/functions_root_pol.m";

/* #################################################### */
/* #################################################### */
/* #################################################### */

SetClassGroupBounds("GRH");
p<x> := PolynomialRing(Integers());
q<y> := PolynomialRing(Integers());

if TYPE_FIELD eq 0 then
    TYPE_FIELD := "complex";
elif TYPE_FIELD eq 1 then
    TYPE_FIELD := "real";
end if;

PRECISIONS := [50, 100, 500, 1000, 2000, 5000];

for i in [0..NUMBER_DIM-1] do
    
    field_dim := DIM_m + i*5;
    if TYPE_FIELD eq "complex" then
	field_dim := 2*(field_dim div 2);
    end if;

    print "********** DIM IS:", field_dim, " **********";

    prune_list := [RealField()!i/field_dim : i in [1..field_dim]];
    string_dim := Sprintf("../data/gh_%o_%o", field_dim, TYPE_FIELD);

    for precision in PRECISIONS do
	QUO := 0;
	
	printf "field_dim: %3o, PRECISION: %3o \n", field_dim, precision;
	
	for j in [1..NUMBER_TESTS] do
	    
	    if (j mod 20) eq 0 then
		print j;
	    end if;
	    precision_field := precision + 3*field_dim;
	    
	    size_pol_field := Random(1,10);
	    if TYPE_FIELD eq "real" then
		pol_field := Pol_field_creation_real(field_dim, size_pol_field);
	    elif TYPE_FIELD eq "complex" then
		pol_field := Pol_field_creation_complex(field_dim, size_pol_field);
	    end if;
	    
	    K := NumberField(pol_field: Check:=false);
	    a := K.1;
	    r1, r2 := Signature(K);
	    
	    if TYPE_FIELD eq "real" then
		B := RealEmbeddings(a);
		B := [Abs(B[i]) : i in [1..#B]];
		m, ind := Max(B);
		B1 := Create_real_matrix_basis(K, precision, precision_field:
					       index:=ind);
	    elif TYPE_FIELD eq "complex" then
		B := Conjugates(a)[r1+1..r1+2*r2];
		B := [Abs(B[i]) : i in [1..#B]];
		m, ind := Max(B);
		ind := r1+ind;
		B1 := Create_complex_matrix_basis(K, precision: index:=ind);
	    end if;
	    
	    L := MyLLL(B1: version:="fplll", Delta:=0.99);
	    rank := Rank(L);
	    vol := Sqrt(Determinant(L*Transpose(L)));
	    
	    l1 := L1_est(vol, rank: est:="gh");
	    
	    sv_n := Sqrt(Norm(ShortestVectors(Lattice(L): Proof:=false, Max:=1, Prune:=prune_list)[1]));

	    QUO +:= sv_n/l1;
	    
	end for;
	print QUO/NUMBER_TESTS;
	FILE := Open(string_dim, "a");
	fprintf FILE, "%o \t %o\n", precision, RealField(6)!QUO/NUMBER_TESTS;
	delete FILE;

    end for;

end for;
