read("../src/functions_root_pol_gp.c");
read("../src/root_pol_multiemb.c");

TYPE_FIELD = Str(TYPE_FIELD);


default(parisizemax, 1000000000000);
{
  my(n, Q, p, quo, str_quo, FILE_PRINT, g, prec, prec_t, latt, uni, base, h,\
     gemb, SIZE_ROOTS_VEC, SIZE_ROOTS);

  SIZE_ROOTS_VEC = [1, 25, 50, 75, 100];

  for (k = 1, length(SIZE_ROOTS_VEC),
	 SIZE_ROOTS = SIZE_ROOTS_VEC[k];


       printf("*************************************************** \n");
       printf("************* SIZE of ROOTS IS =  %s  *********** \n", SIZE_ROOTS);
       FILE_PRINT = strprintf("../data/prec_eval_%s_%s_%s", TYPE_FIELD, SIZE_POL_FIELD, SIZE_ROOTS);
  
       for (i = 0, NUMBER_DIM-1,
	      n = DIM_START + i*25;
	    if (TYPE_FIELD=="complex", n = 2*(n \ 2));
	    printf("************* DIM =  %s  *********** \n", n);
       
	    Q = [];
	    for(j = 1, NUMBER_TESTS,
		  if (j%20==0,
		      printf("        *** %s th test *** \n", j);
		      );
		if (TYPE_FIELD=="real",
		    p = Pol_field_creation_real(n, SIZE_POL_FIELD); ,
		    p = Pol_field_creation_complex(n, SIZE_POL_FIELD);
		    );
		g = vector(n, k, random(2*2^SIZE_ROOTS)-2^SIZE_ROOTS);
		prec_t = n*log(norml2(g))/(2*log(2));
	   
		h = 0;
		c = 0;
		prec = max(round(prec_t)-200, 10);
		if (TYPE_FIELD=="real", 
		    [latt, uni, base, prec] = New_basis(p, prec, prec, 0, 0); ,
		    [latt, uni, base, prec] = New_basis_complex(p, prec, prec, 0, 0);
		    );

		while(h != g,
		      c +=1;
		      prec += 50;
		      if (TYPE_FIELD=="real", 
			  [latt, uni, base] = Real_basis_lattice(p, prec, uni, \
								 prec, 0); ,
			  [latt, uni, base] = Complex_basis_lattice(p, prec, \
								    uni, prec, 0);
			  );
		      gemb = Element_embedding(g, base);
		      if (TYPE_FIELD=="real",
			  h = Test_decode([gemb], prec, latt); ,
			  TYPE_FIELD=="complex",
			  h = Test_decode([real(gemb), imag(gemb)], prec, latt);
			  );
		      );
		quo = prec/prec_t;
		Q = concat(Q, quo);
		);

	    Q = vecsort(Q);
	    default(realprecision, 7);
	    str_quo = strprintf("%s\t%s\t%s\t%s", n, Q[1], Q[length(Q)],\
				vecsum(Q)/length(Q));
	    f = fileopen(FILE_PRINT, "a"); /* printing for cert method */
	    filewrite(f, str_quo);
	    fileclose(f);
	    );
       );
};
exit;
