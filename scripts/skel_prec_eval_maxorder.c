read("../src/functions_root_pol_gp.c");
read("../src/root_pol_multiemb.c");

TYPE_FIELD = Str(TYPE_FIELD);
/* TYPE_FIELD = "real"; */

/* SIZE_POL_FIELD=1; */
/* SIZE_ROOTS=10; */
/* NUMBER_DIM=2; */
/* DIM_START=15; */
/* NUMBER_TESTS=10; */

default(parisizemax, 1000000000000);
{
  my(n, Q, p, quo, str_quo, FILE_PRINT, g, prec, prec_t, latt, uni, base, h, \
     gemb, SIZE_ROOTS_VEC, SIZE_ROOTS);
  
  SIZE_ROOTS_VEC = [25];
  
  for (kk = 1, length(SIZE_ROOTS_VEC),
	 SIZE_ROOTS = SIZE_ROOTS_VEC[kk];
       
       FILE_PRINT = strprintf("../data/prec_eval_maxorder_%s_%s_%s", TYPE_FIELD, \
			      SIZE_POL_FIELD, SIZE_ROOTS);
       
       for (i = 0, NUMBER_DIM-1,
	      n = DIM_START + i*5;
	    printf("************* DIM =  %s  *********** \n", n);
	    
	    Q = [];
	    for(j = 1, NUMBER_TESTS,
		  if (j%20==0,
		      printf("     *** %s th test *** \n", j);
		      );
		if (TYPE_FIELD=="real",
		    p = Pol_field_creation_real(n, SIZE_POL_FIELD); ,
		    p = Pol_field_creation_complex(n, SIZE_POL_FIELD);
		    );
	   
		my (K = nfinit(p));
		order_basis = K.zk;
	   
		g = vector(n, k, random(2*2^SIZE_ROOTS)-2^SIZE_ROOTS);
		prec_t = n*log(norml2(g))/(2*log(2));
	   
		h = 0;
		c = 0;
		prec = max(round(prec_t)-200, 10);
	   
		if (TYPE_FIELD=="real", 
		    [latt, uni, base, prec] = New_basis(p, prec, prec, 0, 0,\
							order_basis); ,
		    [latt, uni, base, prec] = New_basis_complex(p, prec, prec, 0, \
								0, order_basis);
		    );

		
		while(h != g,
		      c += 1;
		      prec += 50;
		      if (TYPE_FIELD=="real", 
			  [latt, uni, base] = Real_basis_lattice(p, prec, uni, prec, 0,
								 order_basis); ,
			  [latt, uni, base] = Complex_basis_lattice(p, prec, uni, prec, 0,
								    order_basis);
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
	    str_quo = strprintf("%s\t%s\t%s\t%s", n, Q[1], Q[length(Q)], \
				vecsum(Q)/length(Q));
	    f = fileopen(FILE_PRINT, "a"); 
	    filewrite(f, str_quo);
	    fileclose(f);
	    );
       );
};
exit;  
