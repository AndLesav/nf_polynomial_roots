read("../src/functions_root_pol_gp.c");

TYPE_FIELD = Str(TYPE_FIELD);

default(parisizemax, 1000000000000);
{
  my(n, Q, p, eq, N, quo, str_quo, FILE_PRINT, SIZE_ROOTS, SIZE_ROOTS_VEC);
  SIZE_ROOTS_VEC = [1, 10, 25, 50, 100];
  
  for (j = 1, length(SIZE_ROOTS_VEC),
	 SIZE_ROOTS = SIZE_ROOTS_VEC[j];
       FILE_PRINT  = strprintf("../data/norm_eval_%s_%s_%s", TYPE_FIELD, SIZE_POL_FIELD, SIZE_ROOTS);
       for (i = 0, NUMBER_DIM-1,
	      n = DIM_START + i*5;
	    C = max((n-1)*(n+2)\4, 1);
	    printf("************* DIM =  %s  *********** \n", n);

	    Q = [];
	    for(j = 1, NUMBER_TESTS,
		  if (j%20==0,
		      print("************************ \n");
		      print(j);
		      );
		if (TYPE_FIELD="real",
		    p = Pol_field_creation_real(n, SIZE_POL_FIELD); ,
		    p = Pol_field_creation_complex(n, SIZE_POL_FIELD);
		    );
	   
		[eq, N] = Equation_creation_nf_split(p, DEGREE_EQ, SIZE_ROOTS, "uniform", 1);
		n_heur = Norm_eval(p, eq, 0);
	   
		quo = log(vecmax(N))-log(n_heur);
		quo /= 2;
	   
		Q = concat(Q, quo);
		);
	    Q = vecsort(Q);
	    default(realprecision, 7);
	    str_quo = strprintf("%s\t%s\t%s\t%s", n, Q[1], Q[length(Q)],\
				vecsum(Q)/length(Q));
	    f = fileopen(FILE_PRINT, "a"); /* printing for cert method */
	    filewrite(f, str_quo);
	    fileclose(f);
	    printf("************************ \n\n");
	    );
       );
};
exit;
