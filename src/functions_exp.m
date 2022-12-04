/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @                                                                        @ */
/* @                         MAGMA FUNCTIONS FOR EXPERIMENTS                @ */
/* @                                                                        @ */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */


ExperimentsRelDegreeFixed_Magma := procedure(dim_ext, dim_ground_m, dim_ground_M,
					    deg_eq, vector_field:
					     number_tests:=25, max_size_roots:=100,
					     version:="split")
    /* file where data is printed */
    string_data_file := Sprintf("../data/RelDegree_fixed_MAGMA_%o_%o_%o_%o",
				version, dim_ext, dim_ground_M, deg_eq); 

    /* size of roots for experiments */
    SIZES := [1, 10, 25, 50, 75, 100, 200, 500, 1000, 2000]; 
    
    /* define polynomial ring so that precomputed polynomials can be manipulated */
    R<z> := PolynomialRing(Integers());
    
    if version eq "split" then
	nr := deg_eq;
    else
	nr := 1;
    end if;

    for k in [1..#SIZES] do
	size_roots := SIZES[k];

	if size_roots lt max_size_roots+1 then

	    c_good := 0;
	    c_bad := 0;
	    
	    time_mgm_good := 0;
	    time_mgm_bad := 0;

	    for i in [1..number_tests] do
		p := R!vector_field[i];
		K := NumberField(p: Check:=false);
		P<x> := PolynomialRing(K);	
		
		if deg_eq eq 2 or nr eq deg_eq then
		    equation := P!Equation_creation_field_split(K, deg_eq, size_roots);
		else
		    equation := P!Equation_creation_field_split(K, nr, size_roots);
		    equation *:= (P!Equation_creation_field_multiquad(K, (deg_eq-nr) div 2,
								      size_roots));    
		end if;

		b_good := Is_good_field(K: deg_eq:=deg_eq);
			
		t := Cputime();
		S := Roots(equation);
		if b_good then
		    c_good +:= 1;
		    time_mgm_good +:= Cputime(t);
		else
		    c_bad +:= 1;
		    time_mgm_bad +:= Cputime(t);
		end if;
		delete S, equation, K;
		print time_mgm_good;
		print time_mgm_bad;
	    end for;
	    
	    time_mgm_good /:= c_good;
	    if c_bad ne 0 then
		time_mgm_bad /:= c_bad;
	    end if;
	    str_mgm := Sprintf("%o\t%o\t%o\t%o\n", SIZES[k], 1.*time_mgm_good, time_mgm_bad,  c_bad);

	    FILE := Open(string_data_file, "a");
	    fprintf FILE, str_mgm;
	    delete FILE;    
	end if;
    end for;
end procedure;
