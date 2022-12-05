INDEX := Random(156156161);

/* ########################################################### */
/* ######## FUNCTIONS FOR LLL: Fplll + MyLLL ################ */
/* ########################################################## */
Fplll := function(matrix: mat_transf:=false, Delta:=0.99)
    /* mat_name := matrix; */
    random_tail := IntegerToString(Random(10^9)) cat "_" cat IntegerToString(INDEX);
    output_file := "output_fplll" cat random_tail;
    print_file := "print_file" cat random_tail;
    pre_fplll := ".pre_fplll" cat random_tail;
    matrix_name := "mat_name" cat random_tail;
    mat := [ElementToSequence(matrix[i]) : i in [1..Nrows(matrix)]];
    FILE := Open(pre_fplll, "w");
    fprintf FILE, matrix_name cat " := Matrix(%o); \n", mat;
    delete FILE;
    if not mat_transf then
	dim := 0;
    else
	dim := Nrows(matrix);
    end if;
        
    string := Sprintf("../src/script_magma_fplll.sh %o %o %o %o %o %o", pre_fplll, matrix_name, output_file, print_file, dim,  Delta);
    System(string);

    a := Read(print_file);
    mat_name := eval a;
    uni := 0;
    if mat_transf then
	b := Read(print_file cat "_uni");
	uni := eval b;
    end if;
    System("rm -f" cat " " cat pre_fplll cat " " cat output_file cat " " cat print_file cat " " cat print_file cat "_uni");
    return mat_name, uni;
end function;


/* LLL with options: if simple then computes unitary matrix 
version is either magma either fplll */
MyLLL := function(M: mat_transf:=false, Delta := 0.99, version:="fplll")
    if (version eq "magma") then
	if not mat_transf then
	    B := LLL(M: Delta:=Delta, Proof:=false);
	    return B, 0;
	else
	    B, U := LLL(M: Delta:=Delta);
	    return B, U;
	end if;
    elif (version eq "fplll") then
	return Fplll(M: mat_transf:=mat_transf, Delta:=Delta);
    else
	print "wrong version given";
	return false;
    end if;
end function;




/* LLL rec with precomp */
MyLLL_precomp := function(A: Delta:=0.75, version:="fplll", index:=2)
    dim := Nrows(A);
    r := Ncols(A)-dim;
    B := RowSubmatrixRange(A, 1, index);
    C := Max(Round((dim-1)*(dim+2)/4), 1);
    if index gt 1 then
	B[index][index+r+1] := C*Ceiling((Norm(B[index-1])));
    end if;
    B := MyLLL(B: version:=version, Delta:=0.99);
    for i:=index+1 to dim do
	B[i-1][i-1+r+1] := 0;
	B := Matrix(Append(B[1..i-1], A[i]));
	t := Cputime();
	for j:=1 to i-2 do
	    B[i] +:= (B[i-1][j+r]/C)*A[j+1];
	end for;
	B[i] -:= Round((B[i][1]*A[1][1] + B[i][2]*A[1][2])/Norm(A[1]))*A[1];
	if (i lt dim) then
	    B[i][i+r+1] := C*Ceiling(Sqrt(Norm(B[i-1])));
	end if;
	t := Cputime();
	B := MyLLL(B: version:=version, Delta:=0.99);
    end for;
    return(B);
end function;


/* ######################################################### */
/* ######## Creation of a number field basis lattice ######## */
/* ######## use of a real or complex embedding ######## */
/* ######################################################### */

/* real version */
Create_real_matrix_basis := function(K, precision, precision_field: index:=1)
    ZZ := IntegerRing();
    field_dim := Degree(K);
    C := Max(Round((field_dim-1)*(field_dim+2)/4), 1);
    a := K.1;
    B1 := [Conjugate(a^i, index: Precision:=precision_field + precision div 3) :
	   i in [0..field_dim-1]];
    realBasis1 := [Round(2^(precision)*B1[i]) : i in [1..#B1]];
    realBasis1 := -Vector(realBasis1);
    M := Ceiling(C)*IdentityMatrix(Integers(), field_dim);
    M := HorizontalJoin(Transpose(ChangeRing(Matrix(realBasis1), Integers())), M);
    return M;
end function;




/* complex version */
Create_complex_matrix_basis := function(K, precision: index:=1)
    ZZ := IntegerRing();
    field_dim := Degree(K);
    C := Max(Round((field_dim-1)*(field_dim+2)/4), 1);
    a := K.1;
    B1 := [Conjugate(a^i, index: Precision:=precision div 3) : i in [0..field_dim-1]];
    re_B1 := [Round(2^(precision)*Re(B1[i])) : i in [1..#B1]];
    re_B1 := -Vector(re_B1);
    im_B1 := [Round(2^(precision)*Im(B1[i])) : i in [1..#B1]];
    im_B1 := -Vector(im_B1);
    M := C*IdentityMatrix(Integers(), field_dim);
    M := HorizontalJoin(Transpose(ChangeRing(Matrix(im_B1), Integers())), M);
    M := HorizontalJoin(Transpose(ChangeRing(Matrix(re_B1), Integers())), M);
    return M;
end function;



/* ################################################################# */
/* #                     ELEMENTS FUNCTIONS                        # */
/* ################################################################# */

Embedding_element := function(h, K, B, precision)
    return &+ [RealField(precision)!h[i]*B[i] : i in [1..#h]];
end function;


nf_random_elt := function(K, size)
    return K!(&+[Random(-2^size, 2^size)*(K.1)^j : j in [1..Degree(K)]]);
end function;


/* ################################################################# */
/* #                     POLYNOMIAL GENERATIONS                     #*/
/* ################################################################# */


Pol_field_creation_real := function(field_dim, size_coeff)
    p<x> := PolynomialRing(Integers());
    q<y> := PolynomialRing(Integers());
    r1 := 0;
    while r1 eq 0 do
	P := x^field_dim + p![Random(-2^size_coeff, 2^size_coeff) : i in [1..field_dim]];
    	while (not IsIrreducible(P)) do
	    P := x^field_dim + p![Random(-2^size_coeff, 2^size_coeff) : i in [1..field_dim]];
    	end while;
    	K<a> := (AbsoluteField(NumberField(P)));
	P := q!DefiningPolynomial(K);
	r1 := Signature(K);
    end while;
    return P;
end function;


Pol_field_creation_complex := function(field_dim, size_coeff)
    p<x> := PolynomialRing(Integers());
    q<y> := PolynomialRing(Integers());
    r1 := 1;
    while r1 ne 0 do
	P := x^field_dim + p![Random(-2^size_coeff, 2^size_coeff) : i in [1..field_dim]];
	while (not IsIrreducible(P)) do
	    P := x^field_dim + p![Random(-2^size_coeff, 2^size_coeff) : i in [1..field_dim]];
	end while;
	K<a> := (AbsoluteField(NumberField(P)));
	P := q!DefiningPolynomial(K);
	r1 := Signature(K);
    end while;
    return P;
end function;


Pol_field_creation := function(field_dim, size_coeff: supp:=1)
    p<x> := PolynomialRing(Integers());
    q<y> := PolynomialRing(Integers());
    L_P := [];
    S := SetToSequence(RandomSubset(Set([1..field_dim-1]), Round(supp*(field_dim-1))));
    S := S cat [0];
    for j in [1..1] do
	P := x^field_dim + p!&+[Random(-2^size_coeff, 2^size_coeff)*x^i : i in S];
	while (not IsIrreducible(P)) do
	    P := x^field_dim + p!&+[Random(-2^size_coeff, 2^size_coeff)*x^i : i in S];
	end while;
	Include(~L_P, P);
    end for;    
    r1 := 0;
    K<a> := (AbsoluteField(NumberField(L_P: Check:=false)));
    P := q!DefiningPolynomial(K);
    return P, K;
end function;    


/* integral polynomial with degree=deg with random coeff in 
                           [-2^size_coeff, 2^size_coeff]  */
Pol_creation_integral := function(deg, size_coeff)
    p<x> := PolynomialRing(Integers());
    return p![Random(-2^size_coeff, 2^size_coeff) : i in [0..deg]];
end function;


/* polynomial of degree deg with coefficients in K
   coefficients are elements in ZZ[K.1]       */
Pol_creation_field := function(K, deg, size_coeff: supp:=[1,1])
    local p, x;
    p<x> := PolynomialRing(K);
    pol := 0;
    S := RandomSubset(Set([0..deg-1]), Round(supp[1]*deg));
    T := RandomSubset(Set([0..Degree(K)-1]), Round(supp[2]*Degree(K)));
    for i in SetToSequence(S) do
	coeff := K!&+[Random(-2^size_coeff, 2^size_coeff)*K.1^j : j in T ];
	pol +:= coeff*x^i;
    end for;
    return pol;
end function;


Pol_field_creation_tower := function(extensions_dim, size_coeff: supp := [1,1,1],
								 real := true)
    local p, q, x, y;
    p<x> := PolynomialRing(Integers());
    if real then
	L := [* Pol_field_creation(extensions_dim[1], size_coeff: supp:=supp[1]) *];
    else
	L := [* Pol_field_creation_complex(extensions_dim[1], size_coeff) *];
    end if;
    /* print "ground field created"; */
    K := NumberField(L[#L]: Check:=false);
    /* local q, y; */
    q<y> := PolynomialRing(K);
    for i in [2..#extensions_dim] do
	pol := q!Pol_creation_field(K, extensions_dim[i]-1, size_coeff: supp:=supp[2..3]);
	pol +:= y^extensions_dim[i];
	while (not IsIrreducible(q!pol)) do
	    pol := q!Pol_creation_field(K, extensions_dim[i]-1, size_coeff: supp:=supp[2..3]);
	    pol +:= y^extensions_dim[i];
	    /* pol := Pol_field_creation(extensions_dim[i], size_coeff); */
	end while;
	Append(~L, q!pol);
	K := NumberField(L[#L]: Check:=false);
	q<y> := PolynomialRing(K);
    end for;
    return L, K;
end function;


Pol_field_creation_kummer := function(exp, length, range)
    local p, q, x, y;
    p<x> := PolynomialRing(Integers());
    P := PrimesUpTo(range);
    L := [];
    while #L lt length do
	m := Random(P);
	Include(~L, x^exp-m);
    end while;
    K := NumberField(L: Check:=false);
    return L, K;
end function;



/* polynomial of degree deg with coefficients in K
   coefficients are elements in ZZ[K.1]       */
Equation_creation_field_split := function(K, deg, size_coeff)
    local p, x;
    p<x> := PolynomialRing(K);
    equation := 1;
    for i in [1..deg] do
	coeff := K!(&+[Random(-2^size_coeff, 2^size_coeff)*(K.1)^j : j in [1..Degree(K)]]);
	equation *:= (x - coeff);
    end for;
    return equation;
end function;

/* polynomial of degree deg with coefficients in K
   coefficients are elements in ZZ[K.1]       */
Equation_creation_field_multiquad := function(K, deg, size_coeff)
    local p, x;
    p<x> := PolynomialRing(K);
    equation := 1;
    for i in [1..deg] do
	coeff := K!&+[Random(-2^size_coeff, 2^size_coeff)*K.1^j : j in [1..Degree(K)]];
	equation *:= (x^2 - (coeff^2+1));
    end for;
    return equation;
end function;


/* lambda_1 estimation by Gaussian heuristic */
L1_est := function(vol, rank: est:="gh")
    if est eq "gh" then 
	lambda_1 := Sqrt(rank/(2*Pi(RealField())*Exp(1)))*Root(vol, rank);
    else
	lambda_1 := Root(vol*Gamma(rank/2+1)/(Pi(RealField())^(rank/2)), rank);
    end if;
    return lambda_1;
end function;
