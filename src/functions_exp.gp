/* ############################################################### */
/* ###############  FUNCTIONS FOR EXPERIMENTS  ################### */
/* ############################################################### */

ExperimentsPolField(dim_start, size_pol, size_root, deg_eq, number_tests, {type_field="real", type_eq="split"})=
{
   my(p, S, S_cert, S_heur, equation, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, FILE_GP, time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, times_gp, times_lll_heur, times_lll_cert);
   
   FILE_LLL_CERT = strprintf("./OUTPUTS_BABAI/POLFIELD_LLL_cert_%s_%s", type_field, size_pol);
   FILE_LLL_HEUR = strprintf("./OUTPUTS_BABAI/POLFIELD_LLL_heur_%s_%s", type_field, size_pol);
   FILE_GP  = strprintf("./OUTPUTS_BABAI/POLFIELD_GP_%s_%s", type_field, size_pol); 
   
   for (k = 0, 17,
	
	time_lll_heur = 0;
	time_lll_cert = 0;
	time_gp = 0;
	
	times_gp = [];
	times_lll_heur = [];
	times_lll_cert = [];
	
	c_heur = 0;
	c_cert = 0;
	
	field_dim = dim_start + 5*k;
	if (type_field=="complex",
	    field_dim = 2*(field_dim\2);
	   );
	print("FIELD DIM: " field_dim "\n");
	
	for(i = 1, number_tests,
	   
	   if (i%20==0, print(i "th element \n");
	       print(time_gp " " time_lll_cert " " time_lll_heur);
	      );
	   
	   if (type_field=="real",
	       
	       p = Pol_field_creation_real(field_dim, size_pol); ,
	       
	       p = Pol_field_creation_complex(field_dim, size_pol);
	       
	      );
	   
	     /* print(p); */
	   equation = Equation_creation_nf_split(p, deg_eq, size_root);
	   
	     /* certified lll method */
	   my(s=getabstime());
	   S_cert = Solve_equation(p, equation, 1);
	   time_lll_cert += (getabstime()-s)/1000.0;
	   times_lll_cert = concat(times_lll_cert, (getabstime()-s)/1000.0);
	     /* c_cert += length(S); */
	   
	     /* heuristic lll method */
	   my(s=getabstime());
	   S_heur = Solve_equation(p, equation, 0);
	   time_lll_heur += (getabstime()-s)/1000.0;
	   times_lll_heur = concat(times_lll_heur, (getabstime()-s)/1000.0);
	     /* c_heur += length(S); */
	   
	   
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
	     /* print(time_lll_cert); */
	     /* print(time_lll_heur); */
	     /* print(time_gp); */
	   
	);
	
	D_gp = MyVecSort_two_groups(times_gp);
	D_lll_cert = MyVecSort_two_groups(times_lll_cert);
	D_lll_heur = MyVecSort_two_groups(times_lll_heur);
	
	
	default(realprecision, 10);
	
	write(FILE_LLL_CERT, field_dim "\t" time_lll_cert/number_tests "\t" c_cert/number_tests"\t" D_lll_cert[1] "\t" D_lll_cert[2] "\t" 1.*length(D_lll_cert[4])/number_tests);
	write(FILE_LLL_HEUR, field_dim "\t" time_lll_heur/number_tests "\t" c_heur/number_tests "\t" D_lll_heur[1] "\t" D_lll_heur[2] "\t"  1.*length(D_lll_heur[4])/number_tests);
	write(FILE_GP, field_dim "\t" time_gp/number_tests "\t" D_gp[1] "\t" D_gp[2] "\t"  1.*length(D_gp[4])/number_tests );
       );
   return(1);
};


ExperimentsSizeRoots(dim, size_pol, deg_eq, number_tests, {type_field="real", type_eq="split"})=
{
   my(V_size, p, S, S_cert, S_heur, equation, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, FILE_GP, time_lll_cert, time_lll_heur, time_gp, times_lll_cert, times_lll_heur, times_gp, c_heur, c_cert, D_gp, D_lll_cert, D_lll_heur);
   
   field_dim = dim;
   if (type_field=="complex",
       field_dim = 2*(field_dim\2);
      );
   
   print("FIELD DIM: " field_dim "\n");
   
   FILE_LLL_CERT = strprintf("./OUTPUTS_BABAI/ROOTS_LLL_cert_%s_%s", type_field, field_dim);
   FILE_LLL_HEUR = strprintf("./OUTPUTS_BABAI/ROOTS_LLL_heur_%s_%s", type_field, field_dim);
   FILE_GP  = strprintf("./OUTPUTS_BABAI/ROOTS_GP_%s_%s", type_field, field_dim); 
   
    /* FILE_LLL_CERT = concat(concat("./OUTPUTS_BABAI/ROOTS_LLL_CERT_", type_field), field_dim); */
    /* FILE_LLL_HEUR = concat(concat("./OUTPUTS_BABAI/ROOTS_LLL_HEUR_", type_field), field_dim); */
    /* FILE_GP = concat(concat("./OUTPUTS_BABAI/ROOTS_GP_", type_field), field_dim); */
   
   V_size = [1, 10, 20, 30, 50, 75, 100];
   
   for (k = 1, length(V_size),
	print("Size is: "  V_size[k] "\n");
	times_gp = [];
	times_lll_heur = [];
	times_lll_cert = [];
	
	time_lll_heur = 0;
	time_lll_cert = 0;
	time_gp = 0;
	
	c_heur = 0;
	c_cert = 0;
	
	
	for(i = 1, number_tests,
	   
	   if (i%20==0, print(i "th element \n");
	       print(time_gp " " time_lll_cert " " time_lll_heur);
	      );
	   
	   if (type_field=="real",
	       
	       p = Pol_field_creation_real(field_dim, size_pol); ,
	       
	       p = Pol_field_creation_complex(field_dim, size_pol);
	       
	      );
	   
	   equation = Equation_creation_nf_split(p, deg_eq, V_size[k]);
	   
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
	
	default(realprecision, 10);
	write(FILE_LLL_CERT, V_size[k] "\t" time_lll_cert/number_tests "\t" c_cert/number_tests"\t" D_lll_cert[1] "\t" D_lll_cert[2] "\t" 1.*length(D_lll_cert[4])/number_tests);
	write(FILE_LLL_HEUR, V_size[k] "\t" time_lll_heur/number_tests "\t" c_heur/number_tests "\t" D_lll_heur[1] "\t" D_lll_heur[2] "\t"  1.*length(D_lll_heur[4])/number_tests); 
	write(FILE_GP, V_size[k] "\t" time_gp/number_tests "\t" D_gp[1] "\t" D_gp[2] "\t"  1.*length(D_gp[4])/number_tests ); 
       );
   return(1);
};




/* ExperimentsEquation(dim, size_pol, size_roots, frac_roots, number_tests, {type_field="real", type_eq="split"})= */
/*   { */
   
/*    my(nr, V_size, p, S, equation, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, FILE_GP, time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, total_cert, total_heur); */

/*    field_dim = dim; */

/*    if (type_field=="complex", */
/*        field_dim = 2*(field_dim\2); */
/*        ); */
   
/*    print("FIELD DIM: " field_dim "\n"); */

   
/*    FILE_LLL_CERT = strprintf("./OUTPUTS_BABAI/EQUATION_LLL_cert_%s_%s", type_field, field_dim); */
/*    FILE_LLL_HEUR = strprintf("./OUTPUTS_BABAI/EQUATION_LLL_heur_%s_%s", type_field, field_dim); */
/*    FILE_GP  = strprintf("./OUTPUTS_BABAI/EQUATION_GP_%s_%s", type_field, field_dim);  */


/*     /\* FILE_LLL_CERT = strprintf("./OUTPUTS_BABAI/EQUATION_LLL_cert_%s_%s", type_field, frac_roots); *\/ */
/*    /\* FILE_LLL_HEUR = strprintf("./OUTPUTS_BABAI/EQUATION_LLL_heur_%s_%s", type_field, frac_roots); *\/ */
/*    FILE_GP  = strprintf("./OUTPUTS_BABAI/ROOTS_GP_%s_%s_NEW", type_field, field_dim); */

         
/*    /\* FILE_LLL_CERT = concat(concat("./OUTPUTS_BABAI/EQUATION_LLL_CERT_", type_field), frac_roots); *\/ */
/*    /\* FILE_LLL_HEUR = concat(concat("./OUTPUTS_BABAI/EQUATION_LLL_HEUR_", type_field), frac_roots); *\/ */
/*    /\* FILE_GP = concat(concat("./OUTPUTS_BABAI/EQUATION_GP_", type_field), frac_roots); *\/ */
   

/*    V_deg = [4, 10, 25, 50, 75]; */
   
/*    for (k = 1, length(V_deg), */
	  
/* 	time_lll_heur = 0; */
/* 	time_lll_cert = 0; */
/* 	time_gp = 0; */
	
/* 	c_heur = 0; */
/* 	c_cert = 0; */

/* 	total_heur = 0; */
/* 	total_cert = 0; */
	
/* 	nr = ceil(V_deg[k]*frac_roots/100); */
	
/* 	print("Degree: ", V_deg[k] "\n"); */
/* 	print("There are ", nr, " roots \n"); */
	
/* 	for(i = 1, number_tests, */
	      
/* 	      if (i%20==0, print(i "th element \n"); ); */
	      
/* 	      if (type_field=="real", */

/* 		  p = Pol_field_creation_real(field_dim, size_pol); , */

/* 		  p = Pol_field_creation_complex(field_dim, size_pol); */
		  
/* 		  ); */
	    
/* 	    /\* print(p); *\/ */
/* 	    equation = Equation_creation_nf_split(p, nr, 10); */
/* 	    equation *= Pol_eq_creation_field(V_deg[k]-nr, 5*(V_deg[k]-nr), p); */
	    
/* 	    /\* certified lll method *\/ */
/* 	    my(s=getabstime()); */
/* 	    S = Solve_equation(p, equation, 1); */
/* 	    time_lll_cert += (getabstime()-s)/1000.0; */
/* 	    if (length(S)==nr, c_cert += 1); */
/* 	    total_cert += length(S); */

	    
/* 	    /\* heuristic lll method *\/ */
/* 	    my(s=getabstime()); */
/* 	    S = Solve_equation(p, equation, 0); */
/* 	    time_lll_heur += (getabstime()-s)/1000.0; */
/* 	    /\* print(length(S)); *\/ */
/* 	    if (length(S)==nr, c_heur += 1); */
/* 	    total_heur += length(S); */
	    
/* 	    /\* /\\* algebraic method / gp *\\/ *\/ */
/* 	    /\* my(s=getabstime()); *\/ */
/* 	    /\* S = nfroots(p, equation); *\/ */
/* 	    /\* time_gp += (getabstime()-s)/1000.0; *\/ */
/* 	    ); */
	
/* 	default(realprecision, 5); */
/* 	write(FILE_LLL_CERT, V_deg[k] "\t" time_lll_cert/number_tests "\t" nr "\t" 1.0*c_cert/(number_tests) "\t" 1.0*total_cert/(number_tests)); */
/* 	write(FILE_LLL_HEUR, V_deg[k] "\t" time_lll_heur/number_tests "\t" nr "\t" 1.0*c_heur/(number_tests) "\t" 1.0*total_heur/(number_tests) ); */
/* 	/\* write(FILE_GP, V_deg[k] "\t" time_gp/number_tests); *\/ */
/* 	); */
/*    return(1); */
/*   }; */


ExperimentsNumberSolutions(dim, size_pol, size_roots, frac_roots, number_tests, {type_field="real", type_eq="uniform"})=
  {
    my(nr, V_size, p, S, S_cert, S_heur, equation, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, FILE_GP, time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, D_gp, D_lll_cert, D_lll_heur);
    
    field_dim = dim;
    if (type_field=="complex",
	field_dim = 2*(field_dim\2);
	);
   
    print("FIELD DIM: " field_dim "\n");
    
    FILE_LLL_CERT = strprintf("./OUTPUTS_BABAI/NB_SOLUTIONS_LLL_cert_%s_%s_%s", type_field, zfield_dim, frac_roots);
    FILE_LLL_HEUR = strprintf("./OUTPUTS_BABAI/NB_SOLUTIONS_LLL_HEUR_%s_%s_%s", type_field, field_dim, frac_roots);
    FILE_GP  = strprintf("./OUTPUTS_BABAI/NB_SOLUTIONS_GP_%s_%s_%s", type_field, field_dim, frac_roots);
    
    V_deg = [10, 25, 50, 75, 100];
    
    for (k = 1, length(V_deg),
	   
	   
	   times_gp = [];
	 times_lll_heur = [];
	 times_lll_cert = [];
	 
	 time_lll_heur = 0;
	 time_lll_cert = 0;
	 time_gp = 0;
	 
	 c_heur = 0;
	 c_cert = 0;
	 
	 
	 nr = ceil(V_deg[k]*frac_roots/100);
	 
	 print("Degree: ", V_deg[k] "\n");
	 print("There are ", nr, " roots \n");
	 
	 for (i = 1, number_tests,
	      
		if (i%20==0,
		    print(i "th element \n");
		    default(realprecision, 5);
		    print(time_lll_cert);
		    print(time_lll_heur);
		    print(time_gp);
		    /* print(c_heur); */
		    );
	    
	      if (type_field=="real",
		
		  p = Pol_field_creation_real(field_dim, size_pol); ,
		
		  p = Pol_field_creation_complex(field_dim, size_pol);
		
		  );
	     
	      equation = Equation_creation_nf_split(p, nr, size_roots, type_eq);
	      equation *= Equation_creation_nf_multiquad(p, (V_deg[k]-nr)\2, 2*size_roots);
	      
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
	 	 
	 default(realprecision, 10);
	 write(FILE_LLL_CERT, V_deg[k] "\t" time_lll_cert/number_tests "\t" 1.0*c_cert/number_tests"\t" D_lll_cert[1] "\t" D_lll_cert[2] "\t" 1.*length(D_lll_cert[4])/number_tests);
	 write(FILE_LLL_HEUR, V_deg[k] "\t" time_lll_heur/number_tests "\t" 1.0*c_heur/number_tests "\t" D_lll_heur[1] "\t" D_lll_heur[2] "\t"  1.*length(D_lll_heur[4])/number_tests);
	 write(FILE_GP, V_deg[k] "\t" time_gp/number_tests "\t" D_gp[1] "\t" D_gp[2] "\t"  1.*length(D_gp[4])/number_tests );
	 );
    return(1);
  };


/* ################################################################### */
/* ######################### RELATIVE METHOD ######################### */
/* ################################################################### */


/* ############### COMPARISON CERT vs HEUR methods ############### */

/* comparison between relative methiods : deg of equation */
ExperimentssRelativeCompar_deg(dim_vec, size_pol_vec, size_root, number_tests) =
  {
   my(V_deg, p, S, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, total_cert, total_heur, eq);


   FILE_LLL_CERT = strprintf("./OUTPUTS_BABAI/REL_COMPAR_LLL_cert_%s_%s", dim_vec[1], dim_vec[2]);
   FILE_LLL_HEUR = strprintf("./OUTPUTS_BABAI/REL_COMPAR_LLL_heur_%s_%s", dim_vec[1], dim_vec[2]);
   
   V_deg = [5, 10, 25, 35, 50];
   /* V_deg = [50]; */
   
   for (k = 1, length(V_deg),
	  
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
	
	print("Degree: ", V_deg[k] "\n");
		
	for(i = 1, number_tests,

	    pol_vec = Pol_rel_ext_creation(dim_vec, size_pol_vec);
	    
	    eq = Rel_equation_creation_split(pol_vec, V_deg[k], size_root)[1];

	    /* /\* certified lll method - search + precision *\/ */
	    /* my(s=getabstime()); */
	    /* S = Rel_pol_roots(pol_vec, eq, 1, 1); */
 	    /* time_cert += (getabstime()-s)/1000.0; */
	    /* if (length(S)==V_deg[k], c_cert += 1); */
	    /* total_cert += length(S); */

	    /* /\* partially lll method - search only *\/ */
	    /* my(s=getabstime()); */
	    /* S = Rel_pol_roots(pol_vec, eq, 1, 0); */
 	    /* time_cert_search += (getabstime()-s)/1000.0; */
	    /* if (length(S)==V_deg[k], c_cert_search += 1); */
	    /* total_cert_search += length(S); */

	    /* /\* partially certified lll method - precision only *\/ */
	    /* my(s=getabstime()); */
	    /* S = Rel_pol_roots(pol_vec, eq, 0, 1); */
 	    /* time_cert_prec += (getabstime()-s)/1000.0; */
	    /* if (length(S)==V_deg[k], c_cert_prec += 1); */
	    /* total_cert_prec += length(S); */

	    /* heuristic lll method */
	    my(s=getabstime());
	    S = Rel_pol_roots(pol_vec, eq, 0, 0);
	    time_heur += (getabstime()-s)/1000.0;
	    if (length(S)==V_deg[k], c_heur += 1);
	    total_heur += length(S);
	    
	    if (i%1==0,
		print(i "th element over " number_tests "\n");
		/* print(time_cert); */
		/* print(time_cert_search); */
		/* print(time_cert_prec); */
		print(time_heur);
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
	
	
	str_cert = strprintf("%s\t%s\t%s\t%s\t%s", V_deg[k], time_cert, time_cert_search, c_cert, c_cert_search);
	str_heur = strprintf("%s\t%s\t%s\t%s\t%s", V_deg[k], time_heur, time_cert_prec, c_heur, c_cert_prec);
	
	/* f = fileopen(FILE_LLL_CERT, "a"); /\* printing for cert method *\/ */
	/* filewrite(f, str_cert); */
	/* fileclose(f); */

	/* f = fileopen(FILE_LLL_HEUR, "a"); /\* printing for heur method *\/ */
	/* filewrite(f, str_heur); */
	/* fileclose(f); */
 	);
   return(1);
  };	      



/* comparison between relative methods : size of solutions */
ExperimentsRelativeCompar_sizeroots(dim_vec, size_pol_vec, deg_eq, number_tests) =
  {
    my(V_deg, p, S, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, total_cert, total_heur, eq, SIZES);


   FILE_LLL_CERT = strprintf("./OUTPUTS_BABAI/REL_COMPAR_LLL_cert_prec_%s_%s_%s", deg_eq, dim_vec[1], dim_vec[2]);
   FILE_LLL_HEUR = strprintf("./OUTPUTS_BABAI/REL_COMPAR_LLL_heur_prec_%s_%s_%s", deg_eq, dim_vec[1], dim_vec[2]);
   
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
		print(time_cert);
		print(time_cert_search);
		print(time_cert_prec);
		print(time_heur);
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
   return(1);
  };	      



/* ########################################################################### */

/* ################ GLOBAL vs RELATIVE vs GP :  generic case ################  */



/* FIRST : GLOBAL vs REL number of solutions  */
ExperimentsEquation_rel_abs(dim_vec, vector_field_rel, vector_field_abs, deg_eq, size_roots, number_tests)=
  {
    my(nr, FRACS, S, equation, FILE_LLL_REL, FILE_LLL_ABS, time_lll_rel, time_lll_abs, c_abs, c_rel, total_abs, total_heur, pol_vec);
    
    FILE_LLL_REL = strprintf("./OUTPUTS_BABAI/REL_NBSOL_LLL_rel_%s_%s", dim_vec[1], dim_vec[2]);
    FILE_LLL_ABS = strprintf("./OUTPUTS_BABAI/REL_NBSOL_LLL_abs_%s_%s", dim_vec[1], dim_vec[2]);

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
	 print("Degree: ", FRACS[k] "\n");
	 print("There are ", nr, " roots \n");
	
	 for(i = 1, number_tests,

	       printf ("%d th test\n", i);
	       
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
	 
    return(1);
  };



/* FIRST : GLOBAL vs REL number of solutions  */
ExperimentsRelDegree_rel_abs(dim_ground, deg_eq, vector_field_rel, vector_field_abs, size_roots, number_tests)=
  {
    my(nr, FRACS, S, equation, FILE_LLL_REL, FILE_LLL_ABS, time_lll_rel, time_lll_abs, c_abs, c_rel, total_abs, total_heur, pol_vec, D_gp);
    
    FILE_LLL_REL = strprintf("./OUTPUTS_BABAI/REL_RelDeg_LLL_rel_%s_%s", dim_ground, deg_eq);
    FILE_LLL_ABS = strprintf("./OUTPUTS_BABAI/REL_RelDeg_LLL_abs_%s_%s", dim_ground, deg_eq);
    FILE_GP = strprintf("./OUTPUTS_BABAI/REL_RelDeg_gp_%s_%s", dim_ground, deg_eq);
    
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

	       printf ("%d th test\n", i);
	       
	       if (i%10==0,
		   printf("Time for rel. lll: %s\n", time_lll_rel);
		   printf("Time for abs. lll: %s\n", time_lll_abs);
		   printf("Time for gp: %s\n", time_gp);
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

	 /* f = fileopen(FILE_GP, "a"); /\* printing for GP *\/ */
	 /* filewrite(f, str_abs); */
	 /* fileclose(f); */

	 );

	 
    return(1);
  };



/* ******************************************************************************** */
/* ******************************  CYCLOTOMIC FIELDS ****************************** */
/* ******************************************************************************** */

ExperimentsCyclotomics(dim_max, deg_eq, size_roots, number_tests) = {
  my(p, S, S_cert, S_heur, equation, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, FILE_GP, time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, times_gp, times_lll_heur, times_lll_cert);

  FILE_LLL_ABS = strprintf("./OUTPUTS_BABAI/CYCLO_LLL_abs_%s_new", size_roots);
  FILE_LLL_REL = strprintf("./OUTPUTS_BABAI/CYCLO_LLL_rel_%s_new", size_roots);
  FILE_GP  = strprintf("./OUTPUTS_BABAI/CYCLO_GP_%s_new", size_roots); 

  for (cond = 10, 500,		/* potential conductors */
	 print("***********************************");
       print("Conductor is ", cond);
       dim = eulerphi(cond);
       print("Dimension is ", dim);
       if ((20 < dim) && (dim < dim_max) ,
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
	       times_lll_abs = concat(times_lll_abs, (getabstime()-s)/1000.0);
	    	    
	       /* relative lll method */
	       my(s=getabstime());
	       S_rel = CyclotomicsRelativeRoots(cond, p, equation);
	       time_lll_rel += (getabstime()-s)/1000.0;
	       times_lll_rel = concat(times_lll_rel, (getabstime()-s)/1000.0);

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
	       /* print(time_lll_abs); */
	       /* print(time_lll_rel); */
	       /* print(time_gp); */

	       
	       if (i%20==0, print(i "th element");
		   printf("Retrieved? %d \t %d \n", c_rel, c_abs);
		   printf("Times (gp / abs / rel) are: %s \t %s \t %s \n", time_gp, time_lll_abs, time_lll_rel);
		   );

	    
	       );

	   D_gp = MyVecSort_two_groups(times_gp);
	   D_lll_abs = MyVecSort_two_groups(times_lll_abs);
	   D_lll_rel = MyVecSort_two_groups(times_lll_rel);

	   time_gp /= number_tests;
	   time_lll_abs /= number_tests;
	   time_lll_rel /= number_tests;

	 
	   default(realprecision, 6);
	   
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




ExperimentsEquation_kummer(dim_vec, deg_eq, vector_field, size_roots, number_tests)=
  {
   
   my(nr, V_size, S, equation, FILE_LLL_REL, FILE_LLL_HEUR, time_lll_rel, time_lll_heur, c_heur, c_rel, total_rel, total_heur, pol_vec);

   /* field_dim = dim; */
   
   /* if (type_field=="complex", */
   /*     field_dim = 2*(field_dim\2); */
   /*     ); */
   
   /* print("FIELD DIM: " field_dim "\n"); */

   FILE_LLL_REL = strprintf("./OUTPUTS_BABAI/REL_EQUATION_kummer_rel_%s_%s", dim_vec[1], dim_vec[2]);
   /* FILE_LLL_HEUR = strprintf("./OUTPUTS_BABAI/REL_EQUATION_LLL_heur_%s_%s", dim_vec[1], dim_vec[2]); */
   /* FILE_GP  = strprintf("./OUTPUTS_BABAI/ROOTS_GP_%s_%s_NEW", type_field, field_dim);  */
         
   /* FILE_LLL_CERT = concat(concat("./OUTPUTS_BABAI/EQUATION_LLL_CERT_", type_field), frac_roots); */
   /* FILE_LLL_HEUR = concat(concat("./OUTPUTS_BABAI/EQUATION_LLL_HEUR_", type_field), frac_roots); */
   /* FILE_GP = concat(concat("./OUTPUTS_BABAI/EQUATION_GP_", type_field), frac_roots); */

   V_prop = [1, 25, 50 , 75, 100];
   
   for (k = 1, length(V_prop),
	  
	  time_lll_heur = 0;
	time_lll_rel = 0;
	/* time_gp = 0; */
	
	c_heur = 0;
	c_rel = 0;

	total_heur = 0;
	total_rel = 0;
	
	nr = ceil(deg_eq*V_prop[k]/100);
	nr = max(2*(nr\2), 1);

	print("Degree: ", V_prop[k] "\n");
	print("There are ", nr, " roots \n");
	
	for(i = 1, number_tests,
	      
	      if (i%5==0, print(i "th element \n"); );
	    
	    pol_vec = vector_field[i+(k-1)*number_tests];
	    pol_vec[2] = Mod(pol_vec[2], pol_vec[1]);
	    
	    /* print(p); */
	    /* equation = Equation_creation_nf_split(p, deg_eq, 10); */
	    equation = Rel_equation_creation_split(pol_vec, nr, size_roots)[1];

	    equation = equation*Rel_equation_creation_multiquad(pol_vec, (deg_eq-nr)\2, size_roots*3)[1];
	    
	    /* equation *= Pol_eq_creation_field(V_deg[k]-nr, 5*(V_deg[k]-nr), p); */
	    
	    
	    /* certified lll method */
	    my(s=getabstime());

	    S = Rel_pol_roots(pol_vec, equation, 0);

	    time_lll_rel += (getabstime()-s)/1000.0;

	    if (length(S)==nr, c_rel += 1);

	    total_rel += length(S);

	    kill(S);	    
	    );
	
	default(realprecision, 5);
	write(FILE_LLL_REL, V_prop[k] "\t" time_lll_rel/number_tests "\t" nr "\t" 1.0*c_rel/(number_tests) "\t" 1.0*total_rel/(number_tests));
	/* write(FILE_LLL_HEUR, V_deg[k] "\t" time_lll_heur/number_tests "\t" nr "\t" 1.0*c_heur/(number_tests) "\t" 1.0*total_heur/(number_tests) ); */
	/* write(FILE_GP, V_deg[k] "\t" time_gp/number_tests); */
	);
   return(1);
  };



/* ExperimentsRelativeCompar(dim_vec, size_pol_vec, size_root, number_tests) = */
ExperimentsEquation_kummer_abs(dim_vec, deg_eq, vector_field, size_roots, number_tests)=
  {
   
   my(nr, V_prop, S, equation, FILE_LLL_ABS, time_lll_abs, c_abs, total_abs, p);

   /* FILE_LLL_REL = strprintf("./OUTPUTS_BABAI/REL_EQUATION_kummer_abs_%s_%s", dim_vec[1], dim_vec[2]); */
   FILE_LLL_ABS = strprintf("./OUTPUTS_BABAI/REL_EQUATION_kummer_abs_%s_%s", dim_vec[1], dim_vec[2]);
   /* FILE_GP  = strprintf("./OUTPUTS_BABAI/ROOTS_GP_%s_%s_NEW", type_field, field_dim);  */   

   
   V_prop = [1, 25, 50 , 75, 100];
   
   for (k = 1, length(V_prop),
	  
	  time_lll_abs = 0;
	/* time_lll_cert = 0; */
	
	c_abs = 0;
	/* c_rel = 0; */

	total_abs = 0;
	/* total_rel = 0; */
	
	nr = ceil(deg_eq*V_prop[k]/100);
	nr = max(2*(nr\2), 1);

	print("Degree: ", V_prop[k] "\n");
	print("There are ", nr, " roots \n");
	
	for(i = 1, number_tests,
	      
	      if (i%5==0, print(i "th element \n"); );
	    
	    p = vector_field[i+(k-1)*number_tests];
	    p = subst(p, z, y);

	    equation = Equation_creation_nf_split(p, nr, size_roots);
	    equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-nr)\2, size_roots*3);
	    
	    /*  lll heur method */
	    my(s=getabstime());

	    S = Solve_equation(p, equation, 0);

	    time_lll_abs += (getabstime()-s)/1000.0;

	    if (length(S)==nr, c_abs += 1);

	    total_abs += length(S);

	    kill(S);
	    
	    );
	
	default(realprecision, 5);
	/* write(FILE_LLL_REL, V_prop[k] "\t" time_lll_rel/number_tests "\t" nr "\t" 1.0*c_rel/(number_tests) "\t" 1.0*total_rel/(number_tests)); */
	write(FILE_LLL_ABS, V_prop[k] "\t" time_lll_abs/number_tests "\t" nr "\t" 1.0*c_abs/(number_tests) "\t" 1.0*total_abs/(number_tests) );
	/* write(FILE_GP, V_deg[k] "\t" time_gp/number_tests); */
	);
   return(1);
  };


ExperimentsRoots_kummer(dim_vec, deg_eq, vector_field, number_tests)=
  {
   
   my(nr, V_size, S, equation, FILE_LLL_REL, time_lll_rel, c_rel, total_rel, pol_vec);


   FILE_LLL_REL = strprintf("./OUTPUTS_BABAI/REL_ROOTS_kummer_rel_%s_%s", dim_vec[1], dim_vec[2]);
   /* FILE_LLL_HEUR = strprintf("./OUTPUTS_BABAI/REL_EQUATION_LLL_heur_%s_%s", dim_vec[1], dim_vec[2]); */
   /* FILE_GP  = strprintf("./OUTPUTS_BABAI/ROOTS_GP_%s_%s_NEW", type_field, field_dim);  */
         
   V_size = [1, 10, 25, 50, 75];
   
   for (k = 1, length(V_size),
	  
	  time_lll_heur = 0;
	time_lll_rel = 0;
	/* time_gp = 0; */
	
	c_heur = 0;
	c_rel = 0;

	total_heur = 0;
	total_rel = 0;

	if (deg_eq==2, nr=2; ,
	    nr = 1;
	    );
	
	print("Size: ", V_size[k] "\n");
	print("There are ", nr, " roots \n");
	
	for(i = 1, number_tests,
	      
	      if (i%5==0, print(i "th element \n"); );
	    
	    pol_vec = vector_field[i+(k-1)*number_tests];
	    pol_vec[2] = Mod(pol_vec[2], pol_vec[1]);
	    
	    /* print(p); */

	    if (deg_eq==2, 
		equation = Equation_creation_nf_split(p, deg_eq, V_size[k]); ,
		
		equation = Rel_equation_creation_split(pol_vec, nr, V_size[k])[1];
		equation = equation*Rel_equation_creation_multiquad(pol_vec, (deg_eq-nr)\2, V_size[k]*2)[1];
		);
	    
	    /* equation *= Pol_eq_creation_field(V_deg[k]-nr, 5*(V_deg[k]-nr), p); */
	    
	    
	    /* certified lll method */
	    my(s=getabstime());

	    S = Rel_pol_roots(pol_vec, equation, 0);

	    time_lll_rel += (getabstime()-s)/1000.0;

	    if (length(S)==nr, c_rel += 1);

	    total_rel += length(S);

	    kill(S);	    
	    );
	
	default(realprecision, 5);
	write(FILE_LLL_REL, V_size[k] "\t" time_lll_rel/number_tests "\t" nr "\t" 1.0*c_rel/(number_tests) "\t" 1.0*total_rel/(number_tests));
	/* write(FILE_LLL_HEUR, V_deg[k] "\t" time_lll_heur/number_tests "\t" nr "\t" 1.0*c_heur/(number_tests) "\t" 1.0*total_heur/(number_tests) ); */
	/* write(FILE_GP, V_deg[k] "\t" time_gp/number_tests); */
	);
   return(1);
  };



ExperimentsRoots_kummer_abs(dim_vec, deg_eq, vector_field, number_tests)=
  {
   
   my(nr, V_size, S, equation, FILE_LLL_ABS, time_lll_abs, c_abs, total_abs, p);

   /* FILE_LLL_REL = strprintf("./OUTPUTS_BABAI/REL_EQUATION_kummer_abs_%s_%s", dim_vec[1], dim_vec[2]); */
   FILE_LLL_ABS = strprintf("./OUTPUTS_BABAI/REL_ROOTS_kummer_abs_%s_%s", dim_vec[1], dim_vec[2]);
   /* FILE_GP  = strprintf("./OUTPUTS_BABAI/ROOTS_GP_%s_%s_NEW", type_field, field_dim);  */   

   
   V_size = [1, 10, 25 , 50, 75];
   
   for (k = 1, length(V_size),
	  
	  time_lll_abs = 0;
	/* time_lll_cert = 0; */
	
	c_abs = 0;
	/* c_rel = 0; */

	total_abs = 0;
	/* total_rel = 0; */

	
	
	if (deg_eq==2,

	    nr=2; ,

	    nr = 1;
	    );
	
	print("Size root is: ", V_size[k] "\n");
	print("There are ", nr, " roots \n");
	
	for(i = 1, number_tests,
	      
	      if (i%5==0, print(i "th element \n"); );
	    
	    p = vector_field[i+(k-1)*number_tests];
	    p = subst(p, z, y);

	    if (deg_eq==2, 
	    
		equation = Equation_creation_nf_split(p, nr, V_size[k]); ,

		equation = Equation_creation_nf_split(p, nr, V_size[k]);
		equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-nr)\2, size_roots*2);
		);
	    
	    /*  lll heur method */
	    my(s=getabstime());

	    S = Solve_equation(p, equation, 0);

	    time_lll_abs += (getabstime()-s)/1000.0;

	    if (length(S)==nr, c_abs += 1);

	    total_abs += length(S);

	    kill(S);
	    
	    );
	
	default(realprecision, 5);
	/* write(FILE_LLL_REL, V_prop[k] "\t" time_lll_rel/number_tests "\t" nr "\t" 1.0*c_rel/(number_tests) "\t" 1.0*total_rel/(number_tests)); */
	write(FILE_LLL_ABS, V_size[k] "\t" time_lll_abs/number_tests "\t" nr "\t" 1.0*c_abs/(number_tests) "\t" 1.0*total_abs/(number_tests) );
	/* write(FILE_GP, V_deg[k] "\t" time_gp/number_tests); */
	);
   return(1);
  };



ExperimentsRoots_kummer_gp(dim_vec, deg_eq, vector_field, number_tests)=
{
   my(nr, V_size, S, equation, FILE_GP, time_gp, total_gp, p);
   
   /* FILE_LLL_REL = strprintf("./OUTPUTS_BABAI/REL_EQUATION_kummer_abs_%s_%s", dim_vec[1], dim_vec[2]); */
   /* FILE_LLL_ABS = strprintf("./OUTPUTS_BABAI/REL_EQUATION_kummer_abs_%s_%s", dim_vec[1], dim_vec[2]); */
   FILE_GP  = strprintf("./OUTPUTS_BABAI/REL_ROOTS_kummer_gp_%s_%s", dim_vec[1], dim_vec[2]);   
   
   
   V_size = [1, 10, 25, 50, 75];
   
   for (k = 1, length(V_size),
	
	print(deg_eq);
	
	if (deg_eq==2,
	    nr = 2; ,
	    nr = 1;
	   );
	
	print("Size root is: ", V_size[k] "\n");
	print("There are ", nr, " roots \n");
	
	for(i=1, number_tests,
	   
	   if (i%5==0, print(i "th element \n"); );
	   p = vector_field[i+(k-1)*number_tests];
	   p = subst(p, z, y);
	   
	   if(deg_eq==2,
	      
	      print("case deg = 2");
	      equation = Equation_creation_nf_split(p, nr, V_size[k]); ,
	      
	      equation = Equation_creation_nf_split(p, nr, V_size[k]);
	      equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-nr)\2, size_roots*2);
	      
	     );
	   
	    /*  lll heur method */
	   my(s=getabstime());
	   print("before solving equation");
	   S = nfroots(p, equation);
	   
	   time_gp += (getabstime()-s)/1000.0;
	   
	   kill(S);
	   
	);
	
	default(realprecision, 5);
	/* write(FILE_LLL_REL, V_prop[k] "\t" time_lll_rel/number_tests "\t" nr "\t" 1.0*c_rel/(number_tests) "\t" 1.0*total_rel/(number_tests)); */
	/* write(FILE_LLL_ABS, V_prop[k] "\t" time_lll_abs/number_tests "\t" nr "\t" 1.0*c_abs/(number_tests) "\t" 1.0*total_abs/(number_tests) ); */
	write(FILE_GP, V_size[k] "\t" time_gp/number_tests);
	
       );
   
   return(1);
};
