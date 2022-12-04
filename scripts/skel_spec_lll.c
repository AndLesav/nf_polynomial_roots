read("../src/functions_root_pol_gp.c");
read("../src/root_pol_multiemb.c");

TYPE_FIELD = Str(TYPE_FIELD);


default(parisizemax, 1000000000000);
{
  my(n, Q, p, eq, N, r, Det, quo, str_quo, FILE_PRINT);
  PREC = [1000, 2000, 4000, 6000, 10000, 20000];
  
  for (i = 0, NUMBER_DIM-1,
	 n = DIM_START + i*50;
       if (TYPE_FIELD=="complex",
	   n = 2*(n \ 2);
	   );
       
       FILE_PRINT = strprintf("../data/spec_lll_gp_%s_%s_%s", TYPE_FIELD, n, SIZE_POL_FIELD);
       FILE_PRINT_RAND_HNF = strprintf("../data/spec_lll_rand_hnf_%s_%s_%s", TYPE_FIELD, n, SIZE_POL_FIELD);

       printf("************* DIM =  %s  *********** \n", n);
       for (j = 1, length(PREC),
       
	      prec = PREC[j];
	    printf("    *** prec is %s  *** \n", prec);
	    if (TYPE_FIELD=="real",
		p = Pol_field_creation_real(n, SIZE_POL_FIELD);
		M = Create_basis_lattice(p, prec);
		r = 1;
		Det = vecmax(abs(M[1,]));,
		
		p = Pol_field_creation_complex(n, SIZE_POL_FIELD);
		M = Create_basis_lattice_complex(p, prec);
		r = 2;
		Det = max(vecmax(abs(M[1,])), vecmax(abs(M[2,])));
		);

	    print(n);
	    print(FILE_PRINT);
	    /* *************************** */
	    /* computation for basis lattice */
	    /* *************************** */
	    
	    default(realprecision, 7);
	    my(t = getabstime());
	    L  = qflll(M, 3);
	    t1 = getabstime()-t;
	    print("Normal lll: ", strtime(t1));
	    
	    my(t = getabstime());
	    Ls  = MySpecLLL(M, 2);
	    t1_s = getabstime()-t;
	    default(realprecision, 7);
	    print("Spec lll: ", strtime(t1_s));

	    default(realprecision, 7);
	    str_times = strprintf("%s\t%s\t%s", prec, t1/1000., t1_s/1000.);
	    f = fileopen(FILE_PRINT, "a"); 
	    filewrite(f, str_times);
	    fileclose(f);

	    /* *************************** */
	    /* computation for random hnf */
	    /* *************************** */
	    H = HNF_gen(n, Det, max((n-1)*(n+2)\4, 1), r);

	    default(realprecision, 20);
	    my(t = getabstime());
	    L  = qflll(H, 3);
	    t1 = getabstime()-t;
	    print("Normal lll for rand hnf: ", strtime(t1));
	    
	    my(t = getabstime());
	    Ls  = MySpecLLL(H,  2);
	    t1_s = getabstime()-t;
	    default(realprecision, 20);
	    print("Spec lll for rand hnf: ", strtime(t1_s));

	    default(realprecision, 7);
	    str_times = strprintf("%s\t%s\t%s", prec, t1/1000., t1_s/1000.);
	    f = fileopen(FILE_PRINT_RAND_HNF, "a"); 
	    filewrite(f, str_times);
	    fileclose(f);
	    );
       );
};
exit;
