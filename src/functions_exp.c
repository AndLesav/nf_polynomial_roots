/* ############################################################### */
/* ###############  FUNCTIONS FOR EXPERIMENTS  ################### */
/* ############################################################### */


ExperimentsPolField(dim_start, size_pol, size_root, deg_eq, number_tests, \
		    {number_dim=17, type_field="real", type_eq="split"})=  {
  my(p, S, S_cert, S_heur, equation, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, FILE_GP, \
     time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, times_gp, times_lll_heur, \
     times_lll_cert, D_gp, D_lll_cert, D_lll_heur, str_gp, str_heur, str_cert);
  
  FILE_LLL_CERT = strprintf("../data/POLFIELD_LLL_cert_%s_%s", type_field, size_pol);
  FILE_LLL_HEUR = strprintf("../data/POLFIELD_LLL_heur_%s_%s", type_field, size_pol);
  FILE_GP  = strprintf("../data/POLFIELD_GP_%s_%s", type_field, size_pol); 
  
  for (k = 0, number_dim,
	  
	 time_lll_heur = 0;
       time_lll_cert = 0;
       time_gp = 0;

       times_gp = [];
       times_lll_heur = [];
       times_lll_cert = [];
       
       c_heur = 0;
       c_cert = 0;
       
       field_dim = dim_start + 5*k;
       if (type_field!="real",
	   field_dim = 2*(field_dim\2); /* change dim for totally complex fields */
	   );
       print("FIELD DIM: " field_dim "\n");
	
       for(i = 1, number_tests,

	     if (i%20==0,
		 print(i-1 "th element \n");
		 default(realprecision, 5);
		 printf("Time for cert. lll: %s\n", time_lll_cert/(i-1));
		 printf("Time for heur. lll: %s\n", time_lll_heur/(i-1));
		 printf("Time for gp: %s\n", time_gp/(i-1));
		 );
	   
	   if (type_field=="real",
	       p = Pol_field_creation_real(field_dim, size_pol); ,
	       p = Pol_field_creation_complex(field_dim, size_pol);
	       );
	   
	   equation = Equation_creation_nf_split(p, deg_eq, size_root);

	   /* certified lll method */
	   my(s=getabstime());
	   S_cert = Solve_equation(p, equation, 1);
	   time_lll_cert += (getabstime()-s)/1000.0;
	   times_lll_cert = concat(times_lll_cert, (getabstime()-s)/1000.0);
	   
	   /* heuristic lll method */
	   my(s=getabstime());
	   S_heur = Solve_equation(p, equation, 0);
	   time_lll_heur += (getabstime()-s)/1000.0;
	   times_lll_heur = concat(times_lll_heur, (getabstime()-s)/1000.0);	     
	   
	   /* algebraic method / gp */
	   my(s=getabstime());
	   S = nfroots(p, equation);
	   time_gp += (getabstime()-s)/1000.0;
	   times_gp = concat(times_gp, (getabstime()-s)/1000.0);

	   /* verify if all roots are retrieved */
	   if (vecprod(S)==vecprod(S_cert),
	       c_cert += 1;
	       );
	   if (vecprod(S)==vecprod(S_heur),
	       c_heur += 1;
	       );
	   
	   kill(S);
	   kill(S_cert);
	   kill(S_heur);
	    
	   );
       
       D_gp = MyVecSort_two_groups(times_gp);
       D_lll_cert = MyVecSort_two_groups(times_lll_cert);
       D_lll_heur = MyVecSort_two_groups(times_lll_heur);

       time_lll_cert /=  number_tests;
       time_lll_heur /=  number_tests;
       time_gp /= number_tests;
       
       c_cert /=  1.0*number_tests;
       c_heur /=  1.0*number_tests;
       
       
       default(realprecision, 7);

       str_cert = strprintf("%s\t%s\t%s\t%s\t%s\t%s", field_dim, time_lll_cert, c_cert, D_lll_cert[1], D_lll_cert[2], 1.*length(D_lll_cert[4])/number_tests);
       str_heur = strprintf("%s\t%s\t%s\t%s\t%s\t%s", field_dim, time_lll_heur, c_heur, D_lll_heur[1], D_lll_heur[2], 1.*length(D_lll_heur[4])/number_tests);
       str_gp = strprintf("%s\t%s\t%s\t%s\t%s", field_dim, time_gp, D_gp[1], D_gp[2], 1.*length(D_gp[4])/number_tests);

       
       f = fileopen(FILE_LLL_CERT, "a"); /* printing for cert method */
       filewrite(f, str_cert);
       fileclose(f);
       
       f = fileopen(FILE_LLL_HEUR, "a"); /* printing for heur method */
       filewrite(f, str_heur);
       fileclose(f);

       f = fileopen(FILE_GP, "a"); /* printing for heur method */
       filewrite(f, str_gp);
       fileclose(f);
       
       );
};


ExperimentsSizeRoots(dim, size_pol, deg_eq, number_tests,\
		     {type_field="real", type_eq="split"})= {
  my(SIZES, p, S, S_cert, S_heur, equation, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, \
     FILE_GP, time_lll_cert, time_lll_heur, time_gp, times_lll_cert, times_lll_heur,\
     times_gp, c_heur, c_cert, D_gp, D_lll_cert, D_lll_heur, str_gp, str_heur,\
     str_cert);

  /* change dim for totally complex fields */
  field_dim = dim;
  if (type_field=="complex",
      field_dim = 2*(field_dim\2);
      );
    
  FILE_LLL_CERT = strprintf("../data/SIZEROOTS_LLL_cert_%s_%s_%s", type_field, field_dim, size_pol);
  FILE_LLL_HEUR = strprintf("../data/SIZEROOTS_LLL_heur_%s_%s_%s", type_field, field_dim, size_pol);
  FILE_GP  = strprintf("../data/SIZEROOTS_GP_%s_%s_%s", type_field, field_dim, size_pol); 
  
  SIZES = [1, 10, 20, 30, 50, 75, 100];
  
  for (k = 1, length(SIZES),
	 print("Size is: "  SIZES[k] "\n");

       times_gp = [];
       times_lll_heur = [];
       times_lll_cert = [];
	
       time_lll_heur = 0;
       time_lll_cert = 0;
       time_gp = 0;
	
       c_heur = 0;
       c_cert = 0;

	
       for(i = 1, number_tests,

	     
	     if (i%20==0,
		 print(i-1 "th element \n");
		 default(realprecision, 5);
		 printf("Time for cert. lll: %s\n", time_lll_cert/(i-1));
		 printf("Time for heur. lll: %s\n", time_lll_heur/(i-1));
		 printf("Time for gp: %s\n", time_gp/(i-1));
		 );
	   
	   if (type_field=="real",
	       p = Pol_field_creation_real(field_dim, size_pol); ,
	       p = Pol_field_creation_complex(field_dim, size_pol);
	       );
	    
	   equation = Equation_creation_nf_split(p, deg_eq, SIZES[k]);
	    
	   /* certified lll method */
	   my(s=getabstime());
	   S_cert = Solve_equation(p, equation, 1);
	   time_lll_cert += (getabstime()-s)/1000.0;
	   times_lll_cert = concat(times_lll_cert, (getabstime()-s)/1000.0);
     
	     
	   /* heuristic lll method */
	   my(s=getabstime());
	   S_heur = Solve_equation(p, equation, 0);
	   time_lll_heur += (getabstime()-s)/1000.0;
	   times_lll_heur = concat(times_lll_heur, (getabstime()-s)/1000.0);
	    
	   /* algebraic method / gp */
	   my(s=getabstime());
	   S = nfroots(p, equation);
	   time_gp += (getabstime()-s)/1000.0;
	   times_gp = concat(times_gp, (getabstime()-s)/1000.0);

	   if (vecprod(S)==vecprod(S_cert),
	       c_cert += 1;
	       );
	   if (vecprod(S)==vecprod(S_heur),
	       c_heur += 1;
	       );

	   kill(S);
	   kill(S_cert);
	   kill(S_heur);
	   );

       D_gp = MyVecSort_two_groups(times_gp);
       D_lll_cert = MyVecSort_two_groups(times_lll_cert);
       D_lll_heur = MyVecSort_two_groups(times_lll_heur);

       time_lll_cert /=  number_tests;
       time_lll_heur /=  number_tests;
       time_gp /= number_tests;
       
       c_cert /=  1.0*number_tests;
       c_heur /=  1.0*number_tests;

       
       default(realprecision, 7);
       
       str_cert = strprintf("%s\t%s\t%s\t%s\t%s\t%s", SIZES[k], time_lll_cert, c_cert, D_lll_cert[1], D_lll_cert[2], 1.*length(D_lll_cert[4])/number_tests);
       str_heur = strprintf("%s\t%s\t%s\t%s\t%s\t%s", SIZES[k], time_lll_heur, c_heur, D_lll_heur[1], D_lll_heur[2], 1.*length(D_lll_heur[4])/number_tests);
       str_gp = strprintf("%s\t%s\t%s\t%s\t%s", SIZES[k], time_gp, D_gp[1], D_gp[2], 1.*length(D_gp[4])/number_tests);

       f = fileopen(FILE_LLL_CERT, "a"); /* printing for cert method */
       filewrite(f, str_cert);
       fileclose(f);
       
       f = fileopen(FILE_LLL_HEUR, "a"); /* printing for heur method */
       filewrite(f, str_heur);
       fileclose(f);

       f = fileopen(FILE_GP, "a"); /* printing for heur method */
       filewrite(f, str_gp);
       fileclose(f);
       );
};



ExperimentsNumberSolutions(dim, type_field, frac_roots,\
			   {size_pol=1, size_roots=10, number_tests=50,\
			    type_eq=""})=  {
  my(nr, DEGS, p, S, S_cert, S_heur, equation, field_dim,\
     FILE_LLL_CERT, FILE_LLL_HEUR, FILE_GP, time_lll_cert, time_lll_heur, time_gp,\
     times_lll_cert, times_lll_heur, times_gp, c_heur, c_cert, D_gp, D_lll_cert,\
     D_lll_heur, str_gp, str_heur, str_cert);

  /* change dim for totally complex fields */
  field_dim = dim;
  if (type_field=="complex",
      field_dim = 2*(field_dim\2);
      );
   
  print("FIELD DIM: " field_dim "\n");


  FILE_LLL_CERT = strprintf("../data/NB_SOLUTIONS_LLL_cert_%s_%s_%s", type_field, field_dim, frac_roots);
  FILE_LLL_HEUR = strprintf("../data/NB_SOLUTIONS_LLL_heur_%s_%s_%s", type_field, field_dim, frac_roots);
  FILE_GP  = strprintf("../data/NB_SOLUTIONS_GP_%s_%s_%s", type_field, field_dim, frac_roots); 
    
  DEGS = [10, 25, 50, 75, 100];
    
  for (k = 1, length(DEGS),
	   
	   
	 times_gp = [];
       times_lll_heur = [];
       times_lll_cert = [];
	 
       time_lll_heur = 0;
       time_lll_cert = 0;
       time_gp = 0;
	 
       c_heur = 0;
       c_cert = 0;
	 
	 
       nr = ceil(DEGS[k]*frac_roots/100);
       
       print("Degree: ", DEGS[k], "\n");
       print("There are ", nr, " roots \n");
	 
       for (i = 1, number_tests,

	      if (i%20==0,
		  print(i-1 "th element \n");
		  default(realprecision, 5);
		  printf("Time for cert. lll: %s\n", time_lll_cert/(i-1));
		  printf("Time for heur. lll: %s\n", time_lll_heur/(i-1));
		  printf("Time for gp: %s\n", time_gp/(i-1));
		  );

	    if (type_field=="real",
		p = Pol_field_creation_real(field_dim, size_pol); ,
		p = Pol_field_creation_complex(field_dim, size_pol);
		);
	     
	    equation = Equation_creation_nf_split(p, nr, size_roots, type_eq);
	    equation *= Equation_creation_nf_multiquad(p, (DEGS[k]-nr)\2, 2*size_roots);
	      
	    /* certified lll method */
	    my(s=getabstime());
	    S_cert = Solve_equation(p, equation, 1);
	    time_lll_cert += (getabstime()-s)/1000.0;
	    times_lll_cert = concat(times_lll_cert, (getabstime()-s)/1000.0);

	    
	    /* heuristic lll method */
	    my(s=getabstime());
	    S_heur = Solve_equation(p, equation, 0);
	    time_lll_heur += (getabstime()-s)/1000.0;
	    times_lll_heur = concat(times_lll_heur, (getabstime()-s)/1000.0);
	      
	        
	    /* algebraic method / gp */
	    my(s=getabstime());
	    S = nfroots(p, equation);
	    time_gp += (getabstime()-s)/1000.0;
	    times_gp = concat(times_gp, (getabstime()-s)/1000.0);

	     
	    if (vecprod(S)==vecprod(S_cert), c_cert += 1; );
	    if (vecprod(S)==vecprod(S_heur), c_heur += 1; );
	      
	    kill(S);
	    kill(S_cert);
	    kill(S_heur);
	    );
	 
       D_gp = MyVecSort_two_groups(times_gp);
       D_lll_cert = MyVecSort_two_groups(times_lll_cert);
       D_lll_heur = MyVecSort_two_groups(times_lll_heur);

       time_lll_cert /=  number_tests;
       time_lll_heur /=  number_tests;
       time_gp /= number_tests;
	 
       c_cert /=  1.0*number_tests;
       c_heur /=  1.0*number_tests;
	 
	 
       default(realprecision, 7);

       str_cert = strprintf("%s\t%s\t%s\t%s\t%s\t%s", DEGS[k], time_lll_cert, c_cert, D_lll_cert[1], D_lll_cert[2], 1.*length(D_lll_cert[4])/number_tests);
       str_heur = strprintf("%s\t%s\t%s\t%s\t%s\t%s", DEGS[k], time_lll_heur, c_heur, D_lll_heur[1], D_lll_heur[2], 1.*length(D_lll_heur[4])/number_tests);
       str_gp = strprintf("%s\t%s\t%s\t%s\t%s", DEGS[k], time_gp, D_gp[1], D_gp[2], 1.*length(D_gp[4])/number_tests);



       f = fileopen(FILE_LLL_CERT, "a"); /* printing for cert method */
       filewrite(f, str_cert);
       fileclose(f);
       
       f = fileopen(FILE_LLL_HEUR, "a"); /* printing for heur method */
       filewrite(f, str_heur);
       fileclose(f);

       f = fileopen(FILE_GP, "a"); /* printing for heur method */
       filewrite(f, str_gp);
       fileclose(f);

	 
       );

};


/* ################################################################### */
/* ######################### RELATIVE METHOD ######################### */
/* ################################################################### */


/* ############### COMPARISON CERT vs HEUR methods ############### */

/* comparison between relative methiods : deg of equation */
ExperimentsRelativeCompar_deg(dim_vec, size_pol_vec, size_root, number_tests) = {
  my(DEGS, p, S, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, total_cert, total_heur, eq, str_cert, str_heur);

  FILE_LLL_CERT = strprintf("../data/REL_COMPAR_LLL_degeq_cert_%s_%s_%s", dim_vec[1], dim_vec[2], size_root);
  FILE_LLL_HEUR = strprintf("../data/REL_COMPAR_LLL_degeq_heur_%s_%s_%s", dim_vec[1], dim_vec[2], size_root);
  
  DEGS = [5, 10, 25, 35, 50];
   
  for (k = 1, length(DEGS),
	  
	 time_heur = 0.;	/* prec+search heuristic */
       time_cert = 0.;	/* everything certified */
       time_cert_prec = 0.;	/* only prec certified */
       time_cert_search = 0; /* only search certified */
	
       c_heur = 0;
       c_cert = 0;
       c_cert_prec = 0;
       c_cert_search = 0;

       total_heur = 0;
       total_cert = 0;
       total_cert_prec = 0;
       total_cert_search = 0;
	
       print("Degree: ", DEGS[k] "\n");
		
       for(i = 1, number_tests,

	     pol_vec = Pol_rel_ext_creation(dim_vec, size_pol_vec);
	    
	   eq = Rel_equation_creation_split(pol_vec, DEGS[k], size_root)[1];

	   /* certified lll method - search + precision */
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, eq, 1, 1);
	   time_cert += (getabstime()-s)/1000.0;
	   if (length(S)==DEGS[k], c_cert += 1);
	   total_cert += length(S);

	   /* partially lll method - search only */
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, eq, 1, 0);
	   time_cert_search += (getabstime()-s)/1000.0;
	   if (length(S)==DEGS[k], c_cert_search += 1);
	   total_cert_search += length(S);

	   /* partially certified lll method - precision only */
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, eq, 0, 1);
	   time_cert_prec += (getabstime()-s)/1000.0;
	   if (length(S)==DEGS[k], c_cert_prec += 1);
	   total_cert_prec += length(S);

	   /* heuristic lll method */
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, eq, 0, 0);
	   time_heur += (getabstime()-s)/1000.0;
	   if (length(S)==DEGS[k], c_heur += 1);
	   total_heur += length(S);
	    
	   if (i%20==0,
	       default(realprecision, 10);
	       print(i "th element over " number_tests "\n");
	       print(time_cert/i);
	       print(time_cert_search/i);
	       print(time_cert_prec/i);
	       print(time_heur/i);
	       print("*****************************");
	       );

	   );

       time_cert /=  number_tests;
       time_heur /=  number_tests;
       time_cert_search /=  number_tests;
       time_cert_prec /=  number_tests;
	
       c_cert /=  1.0*number_tests;
       c_heur /=  1.0*number_tests;
       c_cert_search /=  1.0*number_tests;
       c_cert_prec /=  1.0*number_tests;


       default(realprecision, 7);
	
       str_cert = strprintf("%s\t%s\t%s\t%s\t%s", DEGS[k], time_cert, time_cert_search, c_cert, c_cert_search);
       str_heur = strprintf("%s\t%s\t%s\t%s\t%s", DEGS[k], time_heur, time_cert_prec, c_heur, c_cert_prec);
	
       f = fileopen(FILE_LLL_CERT, "a"); /* printing for cert method */
       filewrite(f, str_cert);
       fileclose(f);
       
       f = fileopen(FILE_LLL_HEUR, "a"); /* printing for heur method */
       filewrite(f, str_heur);
       fileclose(f);
       );
};	      



/* comparison between relative methods : size of solutions */
ExperimentsRelativeCompar_sizeroots(dim_vec, size_pol_vec, deg_eq, number_tests) = {
  my(p, S, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, total_cert, total_heur, eq, SIZES, str_cert, str_heur);
  
  FILE_LLL_CERT = strprintf("../data/REL_COMPAR_LLL_sizeroots_cert_%s_%s_%s", deg_eq, dim_vec[1], dim_vec[2]);
  FILE_LLL_HEUR = strprintf("../data/REL_COMPAR_LLL_sizeroots_heur_%s_%s_%s", deg_eq, dim_vec[1], dim_vec[2]);
   
  SIZES = [1, 10, 25, 50, 75, 100];
   
  for (k = 1, length(SIZES),
	 
	 time_heur = 0.;	/* prec+search heuristic */
       time_cert = 0.;	/* everything certified */
       time_cert_prec = 0.;	/* only prec certified */
       time_cert_search = 0; /* only search certified */
	
       c_heur = 0;
       c_cert = 0;
       c_cert_prec = 0;
       c_cert_search = 0;

       total_heur = 0;
       total_cert = 0;
       total_cert_prec = 0;
       total_cert_search = 0;

       printf ("******************************************************* \n");
       printf("*************** Size is: %d ***************\n ", SIZES[k] );
		
       for(i = 1, number_tests,
	     
	     pol_vec = Pol_rel_ext_creation(dim_vec, size_pol_vec);
	   
	   eq = Rel_equation_creation_split(pol_vec, deg_eq, SIZES[k])[1];
	   
	   /* certified lll method - search + precision */
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, eq, 1, 1);
	   time_cert += (getabstime()-s)/1000.0;
	   if (length(S)==deg_eq, c_cert += 1);
	   total_cert += length(S);
	   
	   /* partially lll method - search only */
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, eq, 1, 0);
	   time_cert_search += (getabstime()-s)/1000.0;
	   if (length(S)==deg_eq, c_cert_search += 1);
	   total_cert_search += length(S);
	   
	   /* partially certified lll method - precision only */
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, eq, 0, 1);
	   time_cert_prec += (getabstime()-s)/1000.0;
	   if (length(S)==deg_eq, c_cert_prec += 1);
	   total_cert_prec += length(S);
	   
	   /* heuristic lll method */
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, eq, 0, 0);
	   time_heur += (getabstime()-s)/1000.0;
	   if (length(S)==deg_eq, c_heur += 1);
	   total_heur += length(S);
	   
	   if (i%20==0,
	       print(i "th element over " number_tests "\n");
	       default(realprecision, 5);
	       printf("Time for cert.: %s\n", time_cert/i);
	       printf("Time for cert. search.: %s\n", time_cert_search/i);
	       printf("Time for cert. prec.: %s\n", time_cert_prec/i);
	       printf("Time for heur.: %s\n", time_heur/i);
	       print("*****************************");
	       );

	   );

       time_cert /=  number_tests;
       time_heur /=  number_tests;
       time_cert_search /=  number_tests;
       time_cert_prec /=  number_tests;
	
       c_cert /=  1.0*number_tests;
       c_heur /=  1.0*number_tests;
       c_cert_search /=  1.0*number_tests;
       c_cert_prec /=  1.0*number_tests;


       default(realprecision, 7);
	
       str_cert = strprintf("%s\t%s\t%s\t%s\t%s", SIZES[k], time_cert, time_cert_search, c_cert, c_cert_search);
       str_heur = strprintf("%s\t%s\t%s\t%s\t%s", SIZES[k], time_heur, time_cert_prec, c_heur, c_cert_prec);
	
       f = fileopen(FILE_LLL_CERT, "a"); /* printing for cert method */
       filewrite(f, str_cert);
       fileclose(f);

       f = fileopen(FILE_LLL_HEUR, "a"); /* printing for heur method */
       filewrite(f, str_heur);
       fileclose(f);
       );
};	      



/* ########################################################################### */

/* ################ ABSOLUTE vs RELATIVE :  generic case ################  */


/* FIRST : ABSOLUTE vs RELATIVE number of solutions  */
/* CAREFUL:  NUMBER_TESTS NEEDS TO BE THE SAME AS THE ONE USED TO CREATE FIELDS */
ExperimentsNbSolutions_abs_rel(dim_vec, vector_field_rel, vector_field_abs, deg_eq, size_roots, number_tests)=  {
  my(nr, FRACS, S, equation, FILE_LLL_REL, FILE_LLL_ABS, time_lll_rel, time_lll_abs, c_abs, c_rel, total_abs, total_heur, pol_vec, str_rel, str_abs);
  
  FILE_LLL_REL = strprintf("../data/REL_NBSOL_LLL_rel_%s_%s_%s", dim_vec[1],\
			   dim_vec[2], size_roots);
  FILE_LLL_ABS = strprintf("../data/REL_NBSOL_LLL_abs_%s_%s_%s", dim_vec[1],\
			   dim_vec[2], size_roots);
  
  FRACS = [1, 25, 50, 75, 100];
  
  for (k = 1, length(FRACS),
	 
	 time_lll_abs = 0;
       time_lll_rel = 0;
       
       c_abs = 0;
       c_rel = 0;
       
       total_abs = 0;
       total_rel = 0;
	
       nr = ceil(deg_eq*FRACS[k]/100);
       nr = max(2*(nr\2), 1);

       print("*********************************************\n");
       print("Ratio qf is: ", FRACS[k] "\n");
       print("There are ", nr, " roots \n");
	
       for(i = 1, number_tests,
	       
	   if (i%20==0, print(i "th element \n");
	       print(time_lll_rel);
	       print(time_lll_abs);
	       );

	   /* first relative computation */
	   pol_vec = vector_field_rel[i+(k-1)*number_tests];
	   pol_vec[2] = Mod(pol_vec[2], pol_vec[1]);
	    
	   equation = Rel_equation_creation_split(pol_vec, nr, size_roots)[1];
	   equation = equation*Rel_equation_creation_multiquad(pol_vec, (deg_eq-nr)\2, size_roots*2)[1];
	    	    
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, equation, 0, 0);
	   time_lll_rel += (getabstime()-s)/1000.0;

	   if (length(S)==nr, c_rel += 1);
	   total_rel += length(S);
	     
	   kill(S);
	     
	   /* then absolute computation */
	   p = vector_field_abs[i+(k-1)*number_tests];
	   p = subst(p, z, y);
	     
	   equation = Equation_creation_nf_split(p, nr, size_roots);
	   equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-nr)\2, size_roots*2);
	    
	   my(s=getabstime());
	   S = Solve_equation(p, equation, 0);
	   time_lll_abs += (getabstime()-s)/1000.0;

	   if (length(S)==nr, c_abs += 1);
	   total_abs += length(S);

	   kill(S);

	   );
	

       time_lll_rel /=  number_tests;
       c_rel /=  1.0*number_tests;

       time_lll_abs /=  number_tests;
       c_abs /=  1.0*number_tests;
	 
       default(realprecision, 7);
	 
       str_rel = strprintf("%s\t%s\t%s", FRACS[k], time_lll_rel, c_rel);
       str_abs = strprintf("%s\t%s\t%s", FRACS[k], time_lll_abs, c_abs);

       f = fileopen(FILE_LLL_REL, "a"); /* printing for heur method */
       filewrite(f, str_rel);
       fileclose(f);

       f = fileopen(FILE_LLL_ABS, "a"); /* printing for heur method */
       filewrite(f, str_abs);
       fileclose(f);

       );
};



/* SECOND : ABSOLUTE vs RELATIVE relative degree [L:K]  */
/* CAREFUL:  NUMBER_TESTS NEEDS TO BE THE SAME AS THE ONE USED TO CREATE FIELDS */
ExperimentsRelDegree_abs_rel(dim_ground, deg_eq, vector_field_rel, vector_field_abs, size_roots, number_tests) = {
  my(nr, FRACS, S, equation, FILE_LLL_REL, FILE_LLL_ABS, time_lll_rel, time_lll_abs,\
     c_abs, c_rel, total_abs, total_heur, pol_vec, D_gp, str_rel, str_abs);
  
  FILE_LLL_REL = strprintf("../data/REL_RelDeg_LLL_rel_%s_%s_%s", dim_ground,\
			   deg_eq, size_roots);
  FILE_LLL_ABS = strprintf("../data/REL_RelDeg_LLL_abs_%s_%s_%s", dim_ground,\
			   deg_eq, size_roots);
  
  DIM_EXT = [2, 3, 4, 5];
  
  for (k = 1, length(DIM_EXT),

	 time_lll_rel = 0;
       time_lll_abs = 0;
       time_gp = 0;

       times_gp = [];
	 
       c_abs = 0;
       c_rel = 0;
	 
       total_abs = 0;
       total_rel = 0;
       
       print("*********************************************\n");
       print("Dimension of extension: ", DIM_EXT[k] "\n");
	
       for(i = 1, number_tests,

	       
	   if (i%20==0,
	       printf ("%d th test\n", i-1);
	       printf("Time for rel. lll: %s\n", time_lll_rel/(i-1));
	       printf("Time for abs. lll: %s\n", time_lll_abs/(i-1));
	       );


	   /* first relative computation */
	   pol_vec = vector_field_rel[i+(k-1)*number_tests];
	   pol_vec[2] = Mod(pol_vec[2], pol_vec[1]);
	    
	   equation = Rel_equation_creation_split(pol_vec, deg_eq, size_roots)[1];
	    	    
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, equation, 0, 0);
	   time_lll_rel += (getabstime()-s)/1000.0;

	   if (length(S)==deg_eq, c_rel += 1);
	   /* print(c_rel); */
	   total_rel += length(S);
	     
	   kill(S);
	     
	   /* then absolute computation */
	   p = vector_field_abs[i+(k-1)*number_tests];
	   p = subst(p, z, y);
	   equation = Equation_creation_nf_split(p, deg_eq, size_roots);

	   /* first lll abs. */
	   my(s=getabstime());
	   S = Solve_equation(p, equation, 0);
	   time_lll_abs += (getabstime()-s)/1000.0;

	   if (length(S)==deg_eq, c_abs += 1);
	   /* print(c_abs); */
	   total_abs += length(S);

	   kill(S);

	     
	   /* /\* then gp (abs rep.) *\/ */
	   /* my(s=getabstime()); */
	   /* S = nfroots(p, equation); */
	   /* time_gp += (getabstime()-s)/1000.0; */
	   /* times_gp = concat(times_gp, (getabstime()-s)/1000.0); */

	   /* if (length(S)==nr, c_abs += 1); */
	   /* total_abs += length(S); */

	   /* kill(S); */

	   );
       
       
       time_lll_rel /=  number_tests;
       c_rel /=  1.0*number_tests;
       
       time_lll_abs /=  number_tests;
       c_abs /=  1.0*number_tests;

       /* D_gp = MyVecSort_two_groups(times_gp); */
       
       default(realprecision, 7);
	 
       str_rel = strprintf("%s\t%s\t%s", DIM_EXT[k], time_lll_rel, c_rel);
       str_abs = strprintf("%s\t%s\t%s", DIM_EXT[k], time_lll_abs, c_abs);
       /* str_gp = strprintf("%s\t%s\t%s\t%s", DIM_EXT[k], time_gp, D_gp[1], D_gp[2]); */

	 
       f = fileopen(FILE_LLL_REL, "a"); /* printing for heur method */
       filewrite(f, str_rel);
       fileclose(f);

       f = fileopen(FILE_LLL_ABS, "a"); /* printing for heur method */
       filewrite(f, str_abs);
       fileclose(f);

       );
};



/* ******************************************************************************** */
/* ******************************  CYCLOTOMIC FIELDS ****************************** */
/* ******************************************************************************** */

ExperimentsCyclotomics(dim_max, deg_eq, size_roots, number_tests,\
		       {conds = [20,500]}) = {
  my(p, S, S_cert, S_heur, equation, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, FILE_GP,\
     time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, times_gp, times_lll_heur, times_lll_cert);
  
  FILE_LLL_ABS = strprintf("../data/CYCLO_LLL_abs_%s_%s", deg_eq, size_roots);
  FILE_LLL_REL = strprintf("../data/CYCLO_LLL_rel_%s_%s", deg_eq, size_roots);
  FILE_GP  = strprintf("../data/CYCLO_GP_%s_%s", deg_eq, size_roots); 

  for (cond = conds[1], conds[2],		/* potential conductors */
	 print("***********************************");
       print("Conductor is ", cond);
       dim = eulerphi(cond);
       print("Dimension is ", dim);
       if ((10 < dim) && (dim < dim_max) ,
	   time_lll_abs = 0;
	   time_lll_rel = 0;
	   time_gp = 0;
	   
	   times_gp = [];
	   times_lll_abs = [];
	   times_lll_rel = [];
	   
	   c_abs = 0;
	   c_rel = 0;
	   
	   p = polcyclo(cond);
	   p = subst(p, x, y);
	   
	   for(i = 1, number_tests,
		 equation = Equation_creation_nf_split(p, deg_eq, size_roots);
	       
	       /* absolute lll method */
	       my(s=getabstime());
	       S_abs = Solve_equation(p, equation, 0);
	       time_lll_abs += (getabstime()-s)/1000.0;
	       /* times_lll_abs = concat(times_lll_abs, (getabstime()-s)/1000.0); */
	    	    
	       /* relative lll method */
	       my(s=getabstime());
	       S_rel = CyclotomicsRelativeRoots(cond, p, equation);
	       time_lll_rel += (getabstime()-s)/1000.0;
	       /* times_lll_rel = concat(times_lll_rel, (getabstime()-s)/1000.0); */
	       
	       /* algebraic method / gp */
	       my(s=getabstime());
	       S = nfroots(p, equation);
	       time_gp += (getabstime()-s)/1000.0;
	       times_gp = concat(times_gp, (getabstime()-s)/1000.0);
 
	       if (vecprod(S)==vecprod(S_abs),
		   c_abs += 1;
		   );
	    
	       if (vecprod(S)==vecprod(S_rel),
		   c_rel += 1;
		   );
	       
	       kill(S);
	       kill(S_abs);
	       kill(S_rel);
	       
	       if (i%20==0, print(i "th element");
		   printf("Retrieved? %d \t %d \n", c_rel, c_abs);
		   printf("Times (gp / abs / rel) are: %s \t %s \t %s \n", time_gp/i, time_lll_abs/i, time_lll_rel/i);
		   /* printf("Times (rel) are: %s  \n", time_lll_rel); */
		   );
	       );

	   D_gp = MyVecSort_two_groups(times_gp);
	   /* D_lll_abs = MyVecSort_two_groups(times_lll_abs); */
	   /* D_lll_rel = MyVecSort_two_groups(times_lll_rel); */

	   time_gp /= number_tests;
	   time_lll_abs /= number_tests;
	   time_lll_rel /= number_tests;
	 
	   default(realprecision, 7);
	   
	   str_rel = strprintf("%s\t%s\t%s\t%s", cond, dim, time_lll_rel, c_rel);
	   str_abs = strprintf("%s\t%s\t%s\t%s", cond, dim, time_lll_abs, c_abs);
	   str_gp = strprintf("%s\t%s\t%s\t%s\t%s", cond, dim, time_gp, D_gp[1], D_gp[2]);

	 
	   f = fileopen(FILE_LLL_REL, "a"); /* printing for heur method */
	   filewrite(f, str_rel);
	   fileclose(f);

	   f = fileopen(FILE_LLL_ABS, "a"); /* printing for heur method */
	   filewrite(f, str_abs);
	   fileclose(f);

	   f = fileopen(FILE_GP, "a"); /* printing for heur method */
	   filewrite(f, str_gp);
	   fileclose(f);
	   );
       );
};




/* ******************************************************************************** */
/* ******************************** KUMMER FIELDS ********************************  */
/* ******************************************************************************** */


ExperimentsKummer_LLL_abs(exp, length_seq, deg_eq, vector_field, {number_tests=25, version="split"}) = {
  my(nr, SIZES, S, equation, FILE, p, c_abs, time_lll_abs, size_roots);
  
  FILE = strprintf("../data/KUMMER_LLL_abs_%s_%s_%s", version, exp, length_seq);

  SIZES = [1, 10, 25, 50, 75];

  if (version=="split", nr = deg_eq; , nr = 1;);

  for (k=1, length(SIZES),
	 size_roots = SIZES[k];
       time_lll_abs = 0;
       c_abs=0;
       for(i = 1, number_tests,
	     
	     if (i%20==0, print(i "th element \n"); );
	   
	   p = vector_field[i+(k-1)*number_tests];
	   p = subst(p, z, y);

	   if (deg_eq==2 || nr==deg_eq,
	       equation = Equation_creation_nf_split(p, deg_eq, size_roots); ,

	       equation = Equation_creation_nf_split(p, nr, size_roots);
	       equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-nr)\2, size_roots*2);
	       );
	    
	   my(s=getabstime());
	   S = Solve_equation(p, equation, 0);
	   time_lll_abs += (getabstime()-s)/1000.0;
	   print(time_lll_abs);
	   if (length(S)>=nr, c_abs += 1);
	   /* print(c_abs); */
	   total_abs += length(S);

	   kill(S);
	   if (i%5==0, print(i "th element \n");
	       print(time_lll_abs);
	       );

	   );

       time_lll_abs /=  number_tests;
       c_abs /=  1.0*number_tests;
       
       default(realprecision, 7);
	 
       str_abs = strprintf("%s\t%s\t%s", SIZES[k], time_lll_abs, c_abs);
       
       f = fileopen(FILE, "a"); /* printing for heur method */
       filewrite(f, str_abs);
       fileclose(f);
        
       );
};



ExperimentsKummer_LLL_rel(exp, length_seq, deg_eq, vector_field, number_tests, {version="split"}) = {
  my(nr, SIZES, S, equation, FILE, p, c_rel, time_lll_rel);
    
  FILE = strprintf("../data/KUMMER_LLL_rel_%s_%s_%s", version, exp, length_seq);

  SIZES = [1, 10, 25, 50, 75];
  
  if (version=="split", nr = deg_eq; ,
      nr = 1;
      );
  
  for (k = 1, length(SIZES),
	  
	 time_lll_rel = 0;
       c_rel = 0;
       size_roots = SIZES[k];

       for(i = 1, number_tests,

	     printf ("%d th element\n", i);
	       
	   /* first relative computation */
	   pol_vec = vector_field[i+(k-1)*number_tests];
	   pol_vec[2] = Mod(pol_vec[2], pol_vec[1]);
	   
	   if (deg_eq==2 || nr==deg_eq,
	       equation = Rel_equation_creation_split(pol_vec, deg_eq, size_roots)[1];,
	       equation = Rel_equation_creation_split(pol_vec, nr, size_roots)[1];
	       equation = equation*Rel_equation_creation_multiquad(pol_vec, (deg_eq-nr)\2, size_roots*2)[1];
	       );	   
	   print ("equation created");
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, equation, 0, 0);
	   time_lll_rel += (getabstime()-s)/1000.0;
	   
	   if (length(S)>=nr, c_rel += 1);
	   total_rel += length(S);
	   /* print(c_rel); */
	   kill(S);
	   if (i%5==0, print(i "th element \n");
	       print(time_lll_rel);
	       );
	   );
	

       time_lll_rel /= number_tests;
       c_rel /=  1.0*number_tests;
       print(c_rel);
       default(realprecision, 7);
	 
       str_rel = strprintf("%s\t%s\t%s", SIZES[k], time_lll_rel, c_rel);
       
       f = fileopen(FILE, "a");
       filewrite(f, str_rel);
       fileclose(f);
       
       );
};



ExperimentsKummer_GP(exp, length_seq, deg_eq, vector_field, number_tests, {version="split"}) = {
  my(nr, SIZES, S, equation, FILE, p,  time_gp, times_gp, D_gp);
  
  FILE = strprintf("../data/KUMMER_GP_%s_%s_%s", version, exp, length_seq);

  SIZES = [1, 10, 25, 50, 75];

  if (version=="split", nr = deg_eq; ,
      nr = 1;
      );
  
  
  for (k = 1, length(SIZES),
	 size_roots = SIZES[k];
       
       time_gp = 0;
       times_gp = [];
       for(i = 1, number_tests,
	      
	     if (i%1==0, print(i "th element \n"); );
	   
	   p = vector_field[i+(k-1)*number_tests];
	   p = subst(p, z, y);
	   
	   /* print(p); */
	   /* print(size_roots); */
	   if (deg_eq==2 || nr==deg_eq,
	       equation = Equation_creation_nf_split(p, deg_eq, size_roots); ,
	       
	       equation = Equation_creation_nf_split(p, nr, size_roots);
	       equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-nr)\2, size_roots*2);
	       );	   

	   
	   /* print((equation)); */
	   my(s=getabstime());
	   S = nfroots(p, equation);
	   time_gp += (getabstime()-s)/1000.0;
	   times_gp = concat(times_gp, (getabstime()-s)/1000.0);
	   if (i%1==0,
	       print(time_gp);
	       );
	   kill(S);
	   );
       
       time_gp /=  number_tests;
       D_gp = MyVecSort_two_groups(times_gp);
       default(realprecision, 7);
       
       str_gp = strprintf("%s\t%s\t%s\t%s", SIZES[k], time_gp, D_gp[1], D_gp[2]);
       
       f = fileopen(FILE, "a");
       filewrite(f, str_gp);
       fileclose(f);
       
       );
};

