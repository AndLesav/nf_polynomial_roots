load "../src/functions_root_pol.m";
/* load "../src/general.m"; */

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



for i in [0..NUMBER_DIM-1] do

    field_dim := DIM_START + i*50;
    
    print "############################################################";
    print "         FIELD DIM IS: ", field_dim;
    print "############################################################";
    string_time := Sprintf("../data/spec_lll_%o_%o_%o", TYPE_FIELD, field_dim,
			   SIZE_POL_FIELD);

    for precision in [1000, 2000, 4000, 6000, 10000, 20000] do
	print "        **** prec is ", precision, " ****";

	/* creation of the pol */
	if TYPE_FIELD eq "real" then
	    pol_field := Pol_field_creation_real(field_dim, SIZE_POL_FIELD);
	elif TYPE_FIELD eq "complex" then
	    pol_field := Pol_field_creation_complex(field_dim, SIZE_POL_FIELD);
	end if;
	
	K := NumberField(pol_field);
	a := K.1;
	r1, r2 := Signature(K);
	
	if TYPE_FIELD eq "real" then
	    B := RealEmbeddings(a);
	    B := [Abs(B[i]) : i in [1..#B]];
	    m, ind := Max(B);
	    B_new := Create_real_matrix_basis(K, precision, precision:index:=ind);
	elif TYPE_FIELD eq "complex" then
	    B := Conjugates(a)[r1+1..r1+2*r2];
	    B := [Abs(B[i]) : i in [1..#B]];
	    m, ind := Max(B);
	    ind := r1+ind;
	    B_new := Create_complex_matrix_basis(K, precision: index:=ind);
	end if;

	t := Cputime();
	L1 := MyLLL(B_new: version:="magma", Delta:=0.99);
	tm := Cputime(t);
	print "Time for magma is: ", tm;
	
	t := Cputime();
	L2 := MyLLL(B_new: version:="fplll", Delta:=0.99);
	tf := Cputime(t);
	print "Time for fplll is: ", tf;

	t := Cputime();
	L3 := MyLLL_precomp(B_new: Delta:=0.99, version:="magma");
	tm_s := Cputime(t);
	print "Time for spec magma is: ", tm_s;
	
	t := Cputime();
	L4 := MyLLL_precomp(B_new: Delta:=0.99, version:="fplll");
	tf_s := Cputime(t);
	print "Time for spec fplll is: ", tf_s;

	FILE := Open(string_time, "a");
	fprintf FILE, "%o \t %o \t %o \t %o \t %o \n", precision, RealField(5)!tm, RealField(5)!tm_s, RealField(5)!tf, RealField(5)!tf_s;
	delete FILE;
	
	print "quotient for magma is: ", tm_s/tm;
	print "quotient for fplll is: ", tf_s/tf;
        print "******************************************";

    end for;

end for;



