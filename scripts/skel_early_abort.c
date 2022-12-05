read("../src/functions_root_pol_gp.c");
read("../src/root_pol_multiemb.c");

TYPE_FIELD = Str(TYPE_FIELD);

default(parisizemax, 1000000000000);
{
  my(n, p, NS, nb_steps, str_nbsteps, FILE_PRINT, g, prec, prec_add, latt, uni,
     base, h, gemb, frac_roots, nr, deg, S, SIZE_ROOTS_VEC, SIZE_ROOTS);

  SIZE_ROOTS_VEC = [1, 25, 50, 75, 100];

  for (k = 1, length(SIZE_ROOTS_VEC),
	 SIZE_ROOTS = SIZE_ROOTS_VEC[k];

       printf("*************************************************** \n");
       printf("************* SIZE of ROOTS IS =  %s  *********** \n", SIZE_ROOTS);

       
       FILE_PRINT = strprintf("../data/early_abort_%s_%s_%s", TYPE_FIELD, \
			      SIZE_POL_FIELD, SIZE_ROOTS);
       
       frac_roots = 50;
       deg = 50;
       nr = ceil(deg*frac_roots/100);
       
       for (i = 0, NUMBER_DIM-1,
	      n = DIM_START + i*25;
	    if (TYPE_FIELD=="complex", n = 2*(n\2));
	    printf("************* DIM =  %s  *********** \n", n);
	    
	    NS = [];
	    for(j = 1, NUMBER_TESTS,
		  if (j%20==0,
		      printf("     *** %s th test *** \n", j);
		      );
		if (TYPE_FIELD=="real",
		    p = Pol_field_creation_real(n, SIZE_POL_FIELD); ,
		    p = Pol_field_creation_complex(n, SIZE_POL_FIELD);
		    );

		equation = Equation_creation_nf_split(p, nr, SIZE_ROOTS, "");
		equation *= Equation_creation_nf_multiquad(p, (deg-nr)\2, 2*SIZE_ROOTS);
	   
		[prec, ns, v] = Precision_eval(p, equation, 0);
		prec = round(prec) + ceil(log(log(n))*log(n)*n + (n\2)*log(n)\2);
		prec_add = prec + round(v\2);	   
		if (TYPE_FIELD=="real",
		    [latt, uni, base, prec] = New_basis(p, prec, prec_add);
		    S = Test_solve_equation_real_nbsteps(p,equation,prec,base,latt,\
							 prec_add,0,0);
		    ,
	       
		    [latt, uni, base, prec] = New_basis_complex(p, prec, prec_add);
		    S = Test_solve_equation_complex_nbsteps(p,equation,prec,base,latt,
							    prec_add, 0, 0);
		    );

	   
		if (S[2]!=0, NS = concat(NS, S[2]););
		);

	    NS = vecsort(NS);
	    default(realprecision, 7);
	    str_nbsteps = strprintf("%s\t%s\t%s\t%s", n, NS[1], NS[length(NS)],\
				    1.*vecsum(NS)/length(NS));
	    f = fileopen(FILE_PRINT, "a");
	    filewrite(f, str_nbsteps);
	    fileclose(f);
	    );
       );
};
exit;
