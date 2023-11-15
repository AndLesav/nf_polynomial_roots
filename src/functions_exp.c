/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @                                                                        @ */
/* @                    GP FUNCTIONS FOR EXPERIMENTS                        @ */
/* @                                                                        @ */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */


/* ########################################################################## */
/* #########################  GOOD vs BAD FIELS  ############################ */
/* ########################################################################## */


/* find proportion of good / bad fields for GP nfroots :
 bad field : an inert prime cannot be found in after a few trials  */
ExperimentsProportionBadFields(dim_start, size_pol, number_tests,	\
			       {deg_eq=1, type_field="real", number_dim=15}) = {
  my (FILE, dim, p, c_bad, str);
  FILE = strprintf("../data/proportion_bad_fields_%s_%s", type_field, size_pol);
  dim = dim_start;

  for (i = 0, number_dim - 1,
	 if (type_field == "complex",
	     dim = 2 * (dim \ 2);
	     );
	 print ("***********  dimension is :  ", dim, "  ***********");
	 c_bad = 0;
	 for (j = 1, number_tests,
		
		if (type_field=="real",
		    p = Pol_field_creation_real(dim, size_pol); ,
		    p = Pol_field_creation_complex(dim, size_pol);
		    );

	      if (!(is_good_field(p, deg_eq)[1]), c_bad++);

	      if (j%1000==0, print(j, " th test"));
	      
	      );
       str = strprintf("%s\t%s", dim, 1.*c_bad/number_tests);
       f = fileopen(FILE, "a"); 
       filewrite(f, str);
       fileclose(f);
       print ("*****************************************");
       dim += 5;
       );
};



/* timings for good / bad fields for GP nfroots */
ExperimentsTimingsBadFields(size_pol, deg_eq, {type_field="real", number_tests=20}) = {
  my(p, S, S_cert, S_heur, equation, dim, FILE_LLL_CERT, FILE_LLL_HEUR, FILE_GP, \
     time_lll_cert_good, time_lll_heur_good, time_gp_good, times_gp_bad, \
     times_lll_heur_bad, times_lll_cert_bad, str_gp, str_heur, str_cert);
 
  FILE_LLL_CERT = strprintf("../data/timings_bad_fields_lll_cert_%s_%s", type_field, size_pol);
  FILE_LLL_HEUR = strprintf("../data/timings_bad_fields_lll_heur_%s_%s", type_field, size_pol);
  FILE_GP  = strprintf("../data/timings_bad_fields_gp_%s_%s", type_field, size_pol);

  DIMS = [25, 50, 75, 100];
  
  for (k = 1, 4,
	 
	 dim = DIMS[k];
       /* change dim for totally complex fields */
       if (type_field!="real",  dim = 2*(dim\2); );

       print(" *******    Dimension is : ",  dim, "     *********");
              
       time_lll_heur_good = 0;
       time_lll_cert_good = 0;
       time_gp_good = 0;
       
       time_lll_heur_bad = 0;
       time_lll_cert_bad = 0;
       time_gp_bad = 0;
       
       my(c_good = 0);
       my(c_bad = 0);
       
       while ((c_bad < number_tests) || (c_good < number_tests) ,
	      print(c_good, "   ", c_bad);
	      if (c_good==number_tests,
		  if (type_field=="real",
		      p = Pol_field_creation_real(dim, size_pol); ,
		      p = Pol_field_creation_complex(dim, size_pol);
		      );
		  while (is_good_field(p, deg_eq)[1],
			 if (type_field=="real",
			     p = Pol_field_creation_real(dim, size_pol); ,
			     p = Pol_field_creation_complex(dim, size_pol);
			     );
			 ) ,
		  c_bad==number_tests,
		  if (type_field=="real",
		      p = Pol_field_creation_real(dim, size_pol); ,
		      p = Pol_field_creation_complex(dim, size_pol);
		      );
		  while (!(is_good_field(p, deg_eq)[1]),
			 if (type_field=="real",
			     p = Pol_field_creation_real(dim, size_pol); ,
			     p = Pol_field_creation_complex(dim, size_pol);
			     );
			 ) ,
		  if (type_field=="real",
		      p = Pol_field_creation_real(dim, size_pol); ,
		      p = Pol_field_creation_complex(dim, size_pol);
		      );
		  );

	      b_good = is_good_field(p, deg_eq)[1];
	      if (b_good, c_good++, c_bad++);
	      equation = Equation_creation_nf_split(p, deg_eq, 10);

	      /* certified lll method */
	      my(s=getabstime());
	      S_cert = Solve_equation(p, equation, 1);
	      if (b_good, time_lll_cert_good += (getabstime()-s)/1000.0; ,
		  time_lll_cert_bad += (getabstime()-s)/1000.0; );

	      /* heuristic lll method */
	      my(s=getabstime());
	      S_heur = Solve_equation(p, equation, 0);
	      if (b_good, time_lll_heur_good += (getabstime()-s)/1000.0; ,
		  time_lll_heur_bad += (getabstime()-s)/1000.0; );

	      /* algebraic method / gp */ 
	      my(s=getabstime());
	      S = nfroots(p, equation);
	      if (b_good, time_gp_good += (getabstime()-s)/1000.0; ,
		  time_gp_bad += (getabstime()-s)/1000.0; );
	      /* print(strtime(getabstime() -s )); */
	      
	      kill(S);
	      kill(S_cert);
	      kill(S_heur);
	      
	      );
       
       time_lll_cert_good /= c_good;
       time_lll_heur_good /= c_good;
       time_gp_good /= c_good;
       
       time_lll_cert_bad /= c_bad;
       time_lll_heur_bad /= c_bad;
       time_gp_bad /= c_bad;

       default(realprecision, 7);
  
       /* print results */
       str_cert = strprintf("%s\t%s\t%s", dim, time_lll_cert_good, time_lll_cert_bad);
       str_heur = strprintf("%s\t%s\t%s", dim, time_lll_heur_good, time_lll_heur_bad);
       str_gp = strprintf("%s\t%s\t%s", dim, time_gp_good, time_gp_bad);
      
       
       f = fileopen(FILE_LLL_CERT, "a"); /* printing for cert method */
       filewrite(f, str_cert);
       fileclose(f);
       
       f = fileopen(FILE_LLL_HEUR, "a"); /* printing for heur method */
       filewrite(f, str_heur);
       fileclose(f);

       f = fileopen(FILE_GP, "a"); /* printing for gp method */
       filewrite(f, str_gp);
       fileclose(f);
       
       );
};





/* ########################################################################## */
/* ###########################  ABSOLUTE FIELDS  ############################ */
/* ########################################################################## */



/* experiments for absolute fields : parameters wrt to defining polynomial   */
ExperimentsPolField(dim_start, size_pol, size_root, deg_eq, number_tests, \
		    {number_dim=17, type_field="real", type_eq="split"})=  {

  my(p, S, S_cert, S_heur, equation, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, FILE_GP, \
     time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, times_gp, times_lll_heur, \
     times_lll_cert, D_gp, D_lll_cert, D_lll_heur, str_gp, str_heur, str_cert);
  
  FILE_LLL_CERT = strprintf("../data/POLFIELD_LLL_cert_%s_%s", type_field, size_pol);
  FILE_LLL_HEUR = strprintf("../data/POLFIELD_LLL_heur_%s_%s", type_field, size_pol);
  FILE_GP  = strprintf("../data/POLFIELD_GP_%s_%s", type_field, size_pol);
  
  for (k = 0, number_dim,
	 
	 time_lll_heur_good = 0;
       time_lll_cert_good = 0;
       time_gp_good = 0;
       
       time_lll_heur_bad = 0;
       time_lll_cert_bad = 0;
       time_gp_bad = 0;
       
       c_heur = 0;
       c_cert = 0;

       my(c_good = 0);
       my(c_bad = 0);
       
       field_dim = dim_start + 5*k;
       if (type_field!="real",
	   field_dim = 2*(field_dim\2); /* change dim for totally complex fields */
	   );
       print("FIELD DIM: " field_dim "\n");
	
       for(i = 1, number_tests,

	     if (i%20==0,
		 print(i-1 "th element \n");
	     /* 	 default(realprecision, 5); */
	     /* 	 printf("Time for cert. lll: %s\n", time_lll_cert/(i-1)); */
	     /* 	 printf("Time for heur. lll: %s\n", time_lll_heur/(i-1)); */
	     /* 	 printf("Time for gp: %s\n", time_gp/(i-1)); */
		 );

	   /* randomly generate defining polynomial */
	   if (type_field=="real",
	       p = Pol_field_creation_real(field_dim, size_pol); ,
	       p = Pol_field_creation_complex(field_dim, size_pol);
	       );


	   /* random equation (split) */
	   equation = Equation_creation_nf_split(p, deg_eq, size_root);

	   my(b_good = is_good_field(p, deg_eq)[1]);
	   if (b_good,
	       c_good +=1; ,
	       c_bad += 1;
	       );

	   /* certified lll method */
	   my(s=getabstime());
	   S_cert = Solve_equation(p, equation, 1);
	   if (b_good, time_lll_cert_good += (getabstime()-s)/1000.0; ,
	       time_lll_cert_bad += (getabstime()-s)/1000.0; );
	       	       
	   /* heuristic lll method */
	   my(s=getabstime());
	   S_heur = Solve_equation(p, equation, 0);
	   if (b_good, time_lll_heur_good += (getabstime()-s)/1000.0; ,
	       time_lll_heur_bad += (getabstime()-s)/1000.0; );

	   
	   /* algebraic method / gp */ 
	   my(s=getabstime());
	   S = nfroots(p, equation);
	   if (b_good, time_gp_good += (getabstime()-s)/1000.0; ,
	       time_gp_bad += (getabstime()-s)/1000.0; );

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

       time_lll_cert_good /= c_good;
       time_lll_heur_good /= c_good;
       time_gp_good /= c_good;

       if (c_bad, 
	   time_lll_cert_bad /= c_bad;
	   time_lll_heur_bad /= c_bad;
	   time_gp_bad /= c_bad;
	   );

       
       c_cert /=  1.0*number_tests;
       c_heur /=  1.0*number_tests;

       c_good /=  1.0*number_tests;
       c_bad /=  1.0*number_tests;

       
       default(realprecision, 7);


       /* new print with good vs bad fields */
       str_cert = strprintf("%s\t%s\t%s\t%s\t%s", field_dim, time_lll_cert_good,\
			    time_lll_cert_bad, c_bad, c_cert);
       str_heur = strprintf("%s\t%s\t%s\t%s\t%s", field_dim, time_lll_heur_good, \
			    time_lll_heur_bad, c_bad, c_heur);
       str_gp = strprintf("%s\t%s\t%s\t%s", field_dim, time_gp_good, \
			    time_gp_bad, c_bad);
      
       
       f = fileopen(FILE_LLL_CERT, "a"); /* printing for cert method */
       filewrite(f, str_cert);
       fileclose(f);
       
       f = fileopen(FILE_LLL_HEUR, "a"); /* printing for heur method */
       filewrite(f, str_heur);
       fileclose(f);

       f = fileopen(FILE_GP, "a"); /* printing for gp method */
       filewrite(f, str_gp);
       fileclose(f);
       
       );
};



/* absolute fields : size of roots is varying */
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

       time_lll_heur_good = 0;	/* good and bad : wrt to nfroots (inert prime or not)  */
       time_lll_cert_good = 0;
       time_gp_good = 0;

       time_lll_heur_bad = 0;
       time_lll_cert_bad = 0;
       time_gp_bad = 0;
	
       c_heur = 0; 		/* number of completely solved equations */
       c_cert = 0;

       c_good = 0;		/* number of good / bad fields */
       c_bad = 0;
	
       for(i = 1, number_tests,

	     
	     if (i%20==0,
		 print(i-1 "th element \n");
		 default(realprecision, 5);
		 /* printf("Time for cert. lll: %s\n", time_lll_cert/(i-1)); */
		 /* printf("Time for heur. lll: %s\n", time_lll_heur/(i-1)); */
		 /* printf("Time for gp: %s\n", time_gp/(i-1)); */
		 );

	   /* draw random defining polynomial  */
	   if (type_field=="real",
	       p = Pol_field_creation_real(field_dim, size_pol); ,
	       p = Pol_field_creation_complex(field_dim, size_pol);
	       );

	   /* random equation */
	   equation = Equation_creation_nf_split(p, deg_eq, SIZES[k]);

	   my(b_good = is_good_field(p, deg_eq)[1]);
	   if (b_good, c_good +=1; , c_bad += 1; );
	   
	   /* certified lll method */
	   my(s=getabstime());
	   S_cert = Solve_equation(p, equation, 1);
	   if (b_good, time_lll_cert_good += (getabstime()-s)/1000.0; ,
	       time_lll_cert_bad += (getabstime()-s)/1000.0; );
	     
	   /* heuristic lll method */
	   my(s=getabstime());
	   S_heur = Solve_equation(p, equation, 0);
	   if (b_good, time_lll_heur_good += (getabstime()-s)/1000.0; ,
	       time_lll_heur_bad += (getabstime()-s)/1000.0; );

	    
	   /* algebraic method / gp */
	   my(s=getabstime());
	   S = nfroots(p, equation);
	   if (b_good, time_gp_good += (getabstime()-s)/1000.0; ,
	       time_gp_bad += (getabstime()-s)/1000.0; );

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

       /* D_gp = MyVecSort_two_groups(times_gp); */
       /* D_lll_cert = MyVecSort_two_groups(times_lll_cert); */
       /* D_lll_heur = MyVecSort_two_groups(times_lll_heur); */


       time_lll_cert_good /=  c_good;
       time_lll_heur_good /=  c_good;
       time_gp_good /= c_good;

       if (c_bad, 
	   time_lll_cert_bad /=  c_bad;
	   time_lll_heur_bad /=  c_bad;
	   time_gp_bad /= c_bad;
	   );
       
       c_cert /=  1.0*number_tests;
       c_heur /=  1.0*number_tests;

       c_bad /=  1.0*number_tests;
       c_good /=  1.0*number_tests;
       
       
       
       default(realprecision, 7);
       
       /* str_cert = strprintf("%s\t%s\t%s\t%s\t%s\t%s", SIZES[k], time_lll_cert, c_cert, D_lll_cert[1], D_lll_cert[2], 1.*length(D_lll_cert[4])/number_tests); */
       /* str_heur = strprintf("%s\t%s\t%s\t%s\t%s\t%s", SIZES[k], time_lll_heur, c_heur, D_lll_heur[1], D_lll_heur[2], 1.*length(D_lll_heur[4])/number_tests); */
       /* str_gp = strprintf("%s\t%s\t%s\t%s\t%s", SIZES[k], time_gp, D_gp[1], D_gp[2], 1.*length(D_gp[4])/number_tests); */

       /* new print with good vs bad fields */
       str_cert = strprintf("%s\t%s\t%s\t%s\t%s", SIZES[k], time_lll_cert_good, \
			    time_lll_cert_bad, c_bad, c_cert);
       str_heur = strprintf("%s\t%s\t%s\t%s\t%s", SIZES[k], time_lll_heur_good, \
			    time_lll_heur_bad, c_bad, c_heur);
       str_gp = strprintf("%s\t%s\t%s\t%s", SIZES[k], time_gp_good,	\
			    time_gp_bad, c_bad);
       
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


/* ratio # Z_K(f) / deg f(T)  is varying -- OLD STUFF */
ExperimentsNumberSolutions(dim, type_field, frac_roots,\
			   {size_pol=1, size_roots=10, number_tests=50,\
			    type_eq=""})=  {
  my(nr, DEGS, p, S, S_cert, S_heur, equation, field_dim,		\
     FILE_LLL_CERT, FILE_LLL_HEUR, FILE_GP, time_lll_cert, time_lll_heur, time_gp, \
     times_lll_cert, times_lll_heur, times_gp, c_heur, c_cert, D_gp, D_lll_cert, \
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
	   
	 
	 time_lll_heur_good = 0;
       time_lll_cert_good = 0;
       time_gp_good = 0;
       
       time_lll_heur_bad = 0;
       time_lll_cert_bad = 0;
       time_gp_bad = 0;
	 
       /* times_gp = []; */
       /* times_lll_heur = []; */
       /* times_lll_cert = []; */
	 
       /* time_lll_heur = 0; */
       /* time_lll_cert = 0; */
       /* time_gp = 0; */
       
       c_heur = 0;
       c_cert = 0;

       c_good = 0;
       c_bad = 0;
	 
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
	    equation *= Equation_creation_nf_multiquad(p, (DEGS[k]-nr)\2, size_roots);

	    my(b_good = is_good_field(p, poldegree(equation))[1]);
	    if (b_good, c_good +=1; , c_bad += 1; );
	    
	    /* certified lll method */
	    my(s=getabstime());
	    S_cert = Solve_equation(p, equation, 1);
	    if (b_good, time_lll_cert_good += (getabstime()-s)/1000.0; ,
		time_lll_cert_bad += (getabstime()-s)/1000.0; );

	    /* time_lll_cert += (getabstime()-s)/1000.0; */
	    /* times_lll_cert = concat(times_lll_cert, (getabstime()-s)/1000.0); */

	    
	    /* heuristic lll method */
	    my(s=getabstime());
	    S_heur = Solve_equation(p, equation, 0);
	    if (b_good, time_lll_heur_good += (getabstime()-s)/1000.0; ,
		time_lll_heur_bad += (getabstime()-s)/1000.0; );

	    /* time_lll_heur += (getabstime()-s)/1000.0; */
	    /* times_lll_heur = concat(times_lll_heur, (getabstime()-s)/1000.0); */
	      
	        
	    /* algebraic method / gp */
	    my(s=getabstime());
	    S = nfroots(p, equation);
	    if (b_good, time_gp_good += (getabstime()-s)/1000.0; ,
		time_gp_bad += (getabstime()-s)/1000.0; );

	    /* time_gp += (getabstime()-s)/1000.0; */
	    /* times_gp = concat(times_gp, (getabstime()-s)/1000.0); */

	     
	    if (vecprod(S)==vecprod(S_cert), c_cert += 1; );
	    if (vecprod(S)==vecprod(S_heur), c_heur += 1; );
	      
	    kill(S);
	    kill(S_cert);
	    kill(S_heur);
	    );
	 
       /* D_gp = MyVecSort_two_groups(times_gp); */
       /* D_lll_cert = MyVecSort_two_groups(times_lll_cert); */
       /* D_lll_heur = MyVecSort_two_groups(times_lll_heur); */
       /* time_lll_cert /=  number_tests; */
       /* time_lll_heur /=  number_tests; */
       /* time_gp /= number_tests; */
	 
       time_lll_cert_good /=  c_good;
       time_lll_heur_good /=  c_good;
       time_gp_good /= c_good;

       if (c_bad, 
	   time_lll_cert_bad /=  c_bad;
	   time_lll_heur_bad /=  c_bad;
	   time_gp_bad /= c_bad;
	   );
       
       c_cert /=  1.0*number_tests;
       c_heur /=  1.0*number_tests;

       c_bad /=  1.0*number_tests;
       c_good /=  1.0*number_tests;
       	 
       default(realprecision, 7);


       /* str_cert = strprintf("%s\t%s\t%s\t%s\t%s\t%s", DEGS[k], time_lll_cert, c_cert, D_lll_cert[1], D_lll_cert[2], 1.*length(D_lll_cert[4])/number_tests); */
       /* str_heur = strprintf("%s\t%s\t%s\t%s\t%s\t%s", DEGS[k], time_lll_heur, c_heur, D_lll_heur[1], D_lll_heur[2], 1.*length(D_lll_heur[4])/number_tests); */
       /* str_gp = strprintf("%s\t%s\t%s\t%s\t%s", DEGS[k], time_gp, D_gp[1], D_gp[2], 1.*length(D_gp[4])/number_tests); */

       /* new print with good vs bad fields */
       str_cert = strprintf("%s\t%s\t%s\t%s\t%s", DEGS[k], time_lll_cert_good, \
			    time_lll_cert_bad, c_bad, c_cert);
       str_heur = strprintf("%s\t%s\t%s\t%s\t%s", DEGS[k], time_lll_heur_good, \
			    time_lll_heur_bad, c_bad, c_heur);
       str_gp = strprintf("%s\t%s\t%s\t%s", DEGS[k], time_gp_good,	\
			  time_gp_bad, c_bad);


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




/* ########################################################################## */
/* ########################### RELATIVE METHOD ############################## */
/* ########################################################################## */


/* ##################### COMPARISON CERT vs HEUR methods #################### */

/* comparison between relative methods : deg of equation is varying */
ExperimentsRelativeCompar_deg(dim_vec, size_pol_vec, size_root, number_tests) = {
  my(DEGS, p, S, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, time_lll_cert, \
     time_lll_heur, time_gp, c_heur, c_cert, total_cert, total_heur, eq, \
     str_cert, str_heur);

  FILE_LLL_CERT = strprintf("../data/REL_COMPAR_LLL_degeq_cert_%s_%s_%s", \
			    dim_vec[1], dim_vec[2], size_root);
  FILE_LLL_HEUR = strprintf("../data/REL_COMPAR_LLL_degeq_heur_%s_%s_%s",\
			    dim_vec[1], dim_vec[2], size_root);
  
  DEGS = [5, 10, 25, 35, 50];	/* degree of equations */
   
  for (k = 1, length(DEGS),
	  
	 time_heur = 0.;	/* prec+search heuristic */
       time_cert = 0.;	/* everything certified */
       time_cert_prec = 0.;	/* only prec certified */
       time_cert_search = 0; /* only search certified */
	
       c_heur = 0; /* idem as before : count number of fully solved equations */
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
	
       str_cert = strprintf("%s\t%s\t%s\t%s\t%s", DEGS[k], time_cert, \
			    time_cert_search, c_cert, c_cert_search);
       str_heur = strprintf("%s\t%s\t%s\t%s\t%s", DEGS[k], time_heur, \
			    time_cert_prec, c_heur, c_cert_prec);
	
       f = fileopen(FILE_LLL_CERT, "a"); /* printing for cert method */
       filewrite(f, str_cert);
       fileclose(f);
       
       f = fileopen(FILE_LLL_HEUR, "a"); /* printing for heur method */
       filewrite(f, str_heur);
       fileclose(f);
       );
};	      



/* comparison between relative methods : size of solutions is varying */
ExperimentsRelativeCompar_sizeroots(dim_vec, size_pol_vec, deg_eq, number_tests) = {
  my(p, S, field_dim, FILE_LLL_CERT, FILE_LLL_HEUR, time_lll_cert, time_lll_heur, time_gp, c_heur, c_cert, total_cert, total_heur, eq, SIZES, str_cert, str_heur);
  
  FILE_LLL_CERT = strprintf("../data/REL_COMPAR_LLL_sizeroots_cert_%s_%s_%s",\
			    deg_eq, dim_vec[1], dim_vec[2]);
  FILE_LLL_HEUR = strprintf("../data/REL_COMPAR_LLL_sizeroots_heur_%s_%s_%s",\
			    deg_eq, dim_vec[1], dim_vec[2]);
   
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
	       /* printf("Time for cert.: %s\n", time_cert/i); */
	       /* printf("Time for cert. search.: %s\n", time_cert_search/i); */
	       /* printf("Time for cert. prec.: %s\n", time_cert_prec/i); */
	       /* printf("Time for heur.: %s\n", time_heur/i); */
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
	
       str_cert = strprintf("%s\t%s\t%s\t%s\t%s", SIZES[k], time_cert, \
			    time_cert_search, c_cert, c_cert_search);
       str_heur = strprintf("%s\t%s\t%s\t%s\t%s", SIZES[k], time_heur, \
			    time_cert_prec, c_heur, c_cert_prec);
	
       f = fileopen(FILE_LLL_CERT, "a"); /* printing for cert method */
       filewrite(f, str_cert);
       fileclose(f);

       f = fileopen(FILE_LLL_HEUR, "a"); /* printing for heur method */
       filewrite(f, str_heur);
       fileclose(f);
       );
};



/* ########################################################################## */
/* ################## ABSOLUTE vs RELATIVE :  generic case #################  */
/* ########################################################################## */

/* FIRST : ABSOLUTE vs RELATIVE number of solutions  */
/* CAREFUL:  NUMBER_TESTS NEEDS TO BE THE SAME AS THE ONE USED TO CREATE FIELDS */
ExperimentsNbSolutions_abs_rel(dim_vec, vector_field_rel, vector_field_abs, \
			       deg_eq, size_roots, number_tests)=  {
  my(nr, FRACS, S, equation, FILE_LLL_REL, FILE_LLL_ABS, time_lll_rel, \
     time_lll_abs, c_abs, c_rel, total_abs, total_heur, pol_vec, str_rel, str_abs);
  
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
       /* nr = max(2*(nr\2), 1); */

       print("*********************************************\n");
       print("Ratio qf is: ", FRACS[k] "\n");
       print("There are ", nr, " roots \n");
	
       for(i = 1, number_tests,
	       
	   if (i%20==0, print(i "th element \n");
	       print(time_lll_rel);
	       print(time_lll_abs);
	       );

	   /* first relative computation */
	   pol_vec = vector_field_rel[i];
	   pol_vec[2] = Mod(pol_vec[2], pol_vec[1]);

	   equation = Rel_equation_creation_split(pol_vec, nr, size_roots)[1];
	   /* equation = equation*Rel_pol_creation_mod(pol_vec, deg_eq-nr, size_roots*2); */
	   equation = equation*Rel_equation_creation_multiquad(pol_vec, (deg_eq-nr)\2,size_roots)[1];

	   /* print(equation); */
	   my(s=getabstime());
	   S = Rel_pol_roots(pol_vec, equation, 0, 0);
	   default(realprecision, 12);
	   time_lll_rel += (getabstime()-s)/1000.0;
	   /* printf("time for lll rel is: %s", time_lll_rel); */
	   
	   if (length(S)>=nr, c_rel += 1);
	   total_rel += length(S);
	     
	   kill(S);
	     
	   /* then absolute computation */
	   p = vector_field_abs[i];
	   p = subst(p, z, y);
	   
	   equation = Equation_creation_nf_split(p, nr, size_roots);
	   equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-nr)\2,\
							      size_roots);

	   /* print(equation); */
	   my(s=getabstime());
	   S = Solve_equation(p, equation, 0);
	   time_lll_abs += (getabstime()-s)/1000.0;

	   if (length(S)>=nr, c_abs += 1);
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
/* CAREFUL: NUMBER_TESTS NEEDS TO BE THE SAME AS THE ONE USED TO CREATE FIELDS */
ExperimentsRelDegree_abs_rel(dim_ground, deg_eq, vector_field_rel,\
			     vector_field_abs, size_roots, number_tests) = {
  my(nr, FRACS, S, equation, FILE_LLL_REL, FILE_LLL_ABS, time_lll_rel,\
     time_lll_abs, c_abs, c_rel, total_abs, total_heur, pol_vec, D_gp, str_rel,\
     str_abs);
  
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
	   total_abs += length(S);

	   kill(S);

      	   );
       
       
       time_lll_rel /=  number_tests;
       c_rel /=  1.0*number_tests;
       
       time_lll_abs /=  number_tests;
       c_abs /=  1.0*number_tests;

       
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



/* ########################################################################## */
/* ################## ABSOLUTE vs RELATIVE vs GP vs Magma ################### */
/* ############### One function for each of the possibilities ############### */
/* ##################### relative degree [L:K] fixed ######################## */
/* ########################################################################## */

/* CAREFUL: NUMBER_TESTS NEEDS TO BE SMALLER THAN THE NUMBER OF FIELDS IN */
/*   THE DATABASE FILE USED FOR EXPERIMENTS */
ExperimentsRelDegreeFixed_LLL_abs(dim_ext, dim_ground_m, dim_ground_M, deg_eq, \
				  vector_field, {number_tests=25, max_size_roots=100,\
						 version="split"}) = {
  my(nr, SIZES, S, equation, FILE, p, c_abs, time_lll_abs, size_roots);
  
  FILE = strprintf("../data/RelDegree_fixed_LLL_abs_%s_%s_%s_%s", version, dim_ext,\
		   dim_ground_M, deg_eq);
  
  SIZES = [1, 10, 25, 50, 75, 100, 200, 500, 1000, 2000];

  if (version=="split", nr = deg_eq; , nr = 1;);

  for (k=1, length(SIZES),
	 size_roots = SIZES[k];

       time_lll_abs_good = 0;
       time_lll_abs_bad = 0;

       c_abs=0;

       c_bad=0;
       c_good=0;
       
       if (size_roots < max_size_roots+1,
	   for(i = 1, number_tests,
		 
		 if (i%20==0, print(i "th element \n"); );
	       
	       p = vector_field[2*i-1];
	       p = subst(p, z, y);
	       
	       if (deg_eq==2 || nr==deg_eq,
		   equation = Equation_creation_nf_split(p, deg_eq, size_roots); ,
		   
		   equation = Equation_creation_nf_split(p, nr, size_roots);
		   equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-nr)\2,  size_roots);
		   );

	       my(b_good = is_good_field(p, deg_eq)[1]);
	       if (b_good, c_good +=1; , c_bad += 1; );
	       
	       my(s=getabstime());
	       S = Solve_equation(p, equation, 0);
	       if (b_good, time_lll_abs_good += (getabstime()-s)/1000.0; ,
		   time_lll_abs_bad += (getabstime()-s)/1000.0; );

	       /* time_lll_abs += (getabstime()-s)/1000.0; */
	       /* print(time_lll_abs); */

	       if (length(S)>=nr, c_abs += 1);
	       total_abs += length(S);
	       kill(S);
	       
	       );
	   
	   /* time_lll_abs /=  number_tests; */

	   time_lll_abs_good /=  c_good;

	   if (c_bad, 
	       time_lll_abs_bad /=  c_bad;
	       );

	   c_abs /=  1.0*number_tests;
	   c_bad /=  1.0*number_tests;
	   c_good /=  1.0*number_tests;

	   default(realprecision, 7);
	 
	   str_abs = strprintf("%s\t%s\t%s\t%s\t%s", SIZES[k], time_lll_abs_good, \
			       time_lll_abs_bad, c_bad, c_abs);
       
	   f = fileopen(FILE, "a"); /* printing for heur method */
	   filewrite(f, str_abs);
	   fileclose(f);
	   );
       );
};


ExperimentsRelDegreeFixed_LLL_rel(dim_ext, dim_ground_m, dim_ground_M, deg_eq, \
				  vector_field, {number_tests=25, max_size_roots=100,\
						 version="split"}) = {
  my(nr, SIZES, S, equation, FILE, p, c_rel, time_lll_rel);

  FILE = strprintf("../data/RelDegree_fixed_LLL_rel_%s_%s_%s_%s", version, dim_ext, \
		   dim_ground_M, deg_eq);

  SIZES = [1, 10, 25, 50, 75, 100, 200, 500, 1000, 2000];
  
  if (version=="split", nr = deg_eq; ,
      nr = 1;
      );
  
  for (k = 1, length(SIZES),
	  
	 time_lll_rel_good = 0;
       time_lll_rel_bad = 0;

       c_rel = 0;
       c_bad=0;
       c_good=0;

       size_roots = SIZES[k];

       if (size_roots < max_size_roots+1,
	   for(i = 1, number_tests,

		 printf ("%d th element\n", i);
	   
	       pol_vec = vector_field[i][1..2];
	       b_good = vector_field[i][3];
	       pol_vec[2] = Mod(pol_vec[2], pol_vec[1]);
	       
	       if (deg_eq==2 || nr==deg_eq,
		   equation = Rel_equation_creation_split(pol_vec, deg_eq, size_roots)[1];,
		   equation = Rel_equation_creation_split(pol_vec, nr, size_roots)[1];
		   equation = equation*Rel_equation_creation_multiquad(pol_vec, (deg_eq-nr)\2, size_roots)[1];
		   );

	       
	       /* how to do ?? */
	       /* my(b_good = is_good_field(p, deg_eq)[1]); */
	       if (b_good, c_good +=1; , c_bad += 1; );

	       my(s=getabstime());
	       S = Rel_pol_roots(pol_vec, equation, 0, 0);
	       if (b_good, time_lll_rel_good += (getabstime()-s)/1000.0; ,
		   time_lll_rel_bad += (getabstime()-s)/1000.0; );

	       /* time_lll_rel += (getabstime()-s)/1000.0; */
	   
	       if (length(S)>=nr, c_rel += 1);
	       total_rel += length(S);
	       kill(S);
	       
	       );
       
	   /* time_lll_rel /= number_tests; */

	   time_lll_rel_good /=  c_good;

	   if (c_bad, 
	       time_lll_rel_bad /=  c_bad;
	       );

	   c_rel /=  1.0*number_tests;
	   c_bad /=  1.0*number_tests;
	   c_good /=  1.0*number_tests;

	   default(realprecision, 7);
	 
	   str_rel = strprintf("%s\t%s\t%s\t%s\t%s", SIZES[k], time_lll_rel_good, \
			       time_lll_rel_bad, c_bad, c_rel);
       
	   f = fileopen(FILE, "a");
	   filewrite(f, str_rel);
	   fileclose(f);
	   );       
       );
};

ExperimentsRelDegreeFixed_GP(dim_ext, dim_ground_m, dim_ground_M, deg_eq, \
			     vector_field, {number_tests=25, max_size_roots=100,\
					    version="split"}) = {
  my(nr, SIZES, S, equation, FILE, p,  time_gp, times_gp, D_gp);

  FILE = strprintf("../data/RelDegree_fixed_GP_%s_%s_%s_%s", version, dim_ext, \
		   dim_ground_M, deg_eq);

  SIZES = [1, 10, 25, 50, 75, 100, 200, 500, 1000, 2000];

  if (version=="split", nr = deg_eq; ,
      nr = 1;
      );
  
  
  for (k = 1, length(SIZES),
	 size_roots = SIZES[k];
       
       if (size_roots < max_size_roots+1,
       
	   time_gp_good = 0;
	   time_gp_bad = 0;
	   
	   c_good = 0;
	   c_bad = 0;
	   
	   for(i = 1, number_tests,
	      
		 if (i%1==0, print(i "th element \n"); );
	   
	       p = vector_field[2*i-1];
	       p = subst(p, z, y);
	       
	       if (deg_eq==2 || nr==deg_eq,
		   equation = Equation_creation_nf_split(p, deg_eq, size_roots); ,
	       
		   equation = Equation_creation_nf_split(p, nr, size_roots);
		   equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-nr)\2, size_roots);
		   );	   

	       my(b_good = is_good_field(p, poldegree(equation))[1]);
	       if (b_good, c_good +=1; , c_bad += 1; );

	       
	       /* print((equation)); */
	       my(s=getabstime());
	       S = nfroots(p, equation);
	       if (b_good, time_gp_good += (getabstime()-s)/1000.0; ,
		   time_gp_bad += (getabstime()-s)/1000.0; );

	       if (i%5==0,
		   print(time_gp_good);
		   );
	       kill(S);
	       );

	   time_gp_good /=  c_good;

	   if (c_bad, 
	       time_gp_bad /=  c_bad;
	       );

	   c_bad /=  1.0*number_tests;
	   c_good /=  1.0*number_tests;

	   
	   default(realprecision, 7);
       
	   str_gp = strprintf("%s\t%s\t%s\t%s", SIZES[k], time_gp_good, \
			      time_gp_bad, c_bad);
       
	   f = fileopen(FILE, "a");
	   filewrite(f, str_gp);
	   fileclose(f);
	   );
       );
};









/* ########################################################################## */
/* ########################### CYCLOTOMIC FIELDS ############################ */
/* ########################################################################## */


/* Use the extension Km/Km^+; call on solving function specific to these ext.
 consider all fields with [Km:\Q] = euler_phi(m) < dim_max */
ExperimentsCyclotomics(dim_max, deg_eq, size_roots, number_tests,\
		       {conds = [20,500]}) = {
  my(p, S, S_abs, S_rel, equation, field_dim, FILE_LLL_ABS, FILE_LLL_REL, FILE_GP, \
     time_lll_abs, time_lll_rel, time_gp, c_rel, c_abs, times_gp, times_lll_rel, \
     times_lll_abs);
  
  FILE_LLL_ABS_GOOD = strprintf("../data/CYCLO_LLL_abs_good_%s_%s", deg_eq, size_roots);
  FILE_LLL_REL_GOOD = strprintf("../data/CYCLO_LLL_rel_good_%s_%s", deg_eq, size_roots);
  FILE_GP_GOOD  = strprintf("../data/CYCLO_GP_good_%s_%s", deg_eq, size_roots); 

  FILE_LLL_ABS_BAD = strprintf("../data/CYCLO_LLL_abs_bad_%s_%s", deg_eq, size_roots);
  FILE_LLL_REL_BAD = strprintf("../data/CYCLO_LLL_rel_bad_%s_%s", deg_eq, size_roots);
  FILE_GP_BAD  = strprintf("../data/CYCLO_GP_bad_%s_%s", deg_eq, size_roots); 


  print("deg_eq is: ", deg_eq);
  for (cond = conds[1], conds[2],		/* potential conductors */
	 print("***********************************");
       print("Conductor is ", cond);
       dim = eulerphi(cond);
       print("Dimension is ", dim);
       if ((10 < dim) && (dim < dim_max) ,
	   
	   time_lll_abs = 0;
	   time_lll_rel = 0;
	   time_gp = 0;
	   	   
	   c_abs = 0;
	   c_rel = 0;

	   /* c_good = 0; */
	   /* c_bad = 0; */

	   p = polcyclo(cond);
	   p = subst(p, x, y);
	   print(p);
	   equation = Equation_creation_nf_split(p, deg_eq, size_roots);

	   /* determine if p is a "good" field for this equation size */
	   my(b_good = is_good_field(p, deg_eq)[1]);

	   for(i = 1, number_tests,
		 equation = Equation_creation_nf_split(p, deg_eq, size_roots);
	       
	       /* absolute lll method */
	       my(s=getabstime());
	       S_abs = Solve_equation(p, equation, 0);
	       time_lll_abs += (getabstime()-s)/1000.0;

	       /* relative lll method */
	       my(s=getabstime());
	       S_rel = CyclotomicsRelativeRoots(cond, p, equation);
	       time_lll_rel += (getabstime()-s)/1000.0;

	       /* algebraic method / gp */
	       my(s=getabstime());
	       S = nfroots(p, equation);
	       time_gp += (getabstime()-s)/1000.0;
	       /* times_gp = concat(times_gp, (getabstime()-s)/1000.0); */
	       if (vecprod(S)==vecprod(S_abs),
		   c_abs += 1;
		   );
	    
	       if (vecprod(S)==vecprod(S_rel),
		   c_rel += 1;
		   );
	       
	       kill(S);
	       kill(S_abs);
	       kill(S_rel);
	       
	       if (i%5==0, print(i "th element");
		   printf("Retrieved? %d \t %d \n", c_rel, c_abs);
		   printf("Times (gp / abs / rel) are: %s \t %s \t %s \n", time_gp/i, time_lll_abs/i, time_lll_rel/i);
		   /* printf("Times (rel) are: %s  \n", time_lll_rel); */
		   );
	       );

	   /* D_gp = MyVecSort_two_groups(times_gp); */
	   /* D_lll_abs = MyVecSort_two_groups(times_lll_abs); */
	   /* D_lll_rel = MyVecSort_two_groups(times_lll_rel); */

	   time_gp /= number_tests;
	   time_lll_abs /= number_tests;
	   time_lll_rel /= number_tests;
	 
	   default(realprecision, 7);
	   
	   /* str_rel = strprintf("%s\t%s\t%s\t%s", cond, dim, time_lll_rel, c_rel); */
	   /* str_abs = strprintf("%s\t%s\t%s\t%s", cond, dim, time_lll_abs, c_abs); */
	   /* str_gp = strprintf("%s\t%s\t%s\t%s\t%s", cond, dim, time_gp, D_gp[1], D_gp[2]); */

	   str_rel = strprintf("%s\t%s\t%s\t%s", cond, dim, time_lll_rel, c_rel);
	   str_abs = strprintf("%s\t%s\t%s\t%s", cond, dim, time_lll_abs, c_abs);
	   str_gp = strprintf("%s\t%s\t%s", cond, dim, time_gp);

	   if (b_good, 
	       f = fileopen(FILE_LLL_REL_GOOD, "a"); /* printing for heur method */
	       filewrite(f, str_rel);
	       fileclose(f);
	       
	       f = fileopen(FILE_LLL_ABS_GOOD, "a"); /* printing for heur method */
	       filewrite(f, str_abs);
	       fileclose(f);
	       
	       f = fileopen(FILE_GP_GOOD, "a"); /* printing for heur method */
	       filewrite(f, str_gp);
	       fileclose(f);  ,

	       f = fileopen(FILE_LLL_REL_BAD, "a"); /* printing for heur method */
	       filewrite(f, str_rel);
	       fileclose(f);
	       
	       f = fileopen(FILE_LLL_ABS_BAD, "a"); /* printing for heur method */
	       filewrite(f, str_abs);
	       fileclose(f);
	       
	       f = fileopen(FILE_GP_BAD, "a"); /* printing for heur method */
	       filewrite(f, str_gp);
	       fileclose(f);  
	       
	       );

	   /* f = fileopen(FILE_LLL_REL, "a"); /\* printing for heur method *\/ */
	   /* filewrite(f, str_rel); */
	   /* fileclose(f); */
	   
	   /* f = fileopen(FILE_LLL_ABS, "a"); /\* printing for heur method *\/ */
	   /* filewrite(f, str_abs); */
	   /* fileclose(f); */
	   
	   /* f = fileopen(FILE_GP, "a"); /\* printing for heur method *\/ */
	   /* filewrite(f, str_gp); */
	   /* fileclose(f); */

	   );
       );
};


ExperimentsCyclotomicsNew(dim_max, deg_eq, size_roots, number_tests,	\
			  {conds = [20,500], version="split"}) = {
  my(p, S, S_abs, S_rel, equation, dim, dime, FILE_LLL_ABS, FILE_LLL_REL, FILE_GP, \
     time_lll_abs, time_lll_rel, time_gp, c_rel, c_abs, times_gp, times_lll_rel, \
     times_lll_abs, prime);
  
  FILE_LLL_ABS_GOOD = strprintf("../data/CYCLO_LLL_abs_good_%s_%s_new", deg_eq, size_roots);
  FILE_LLL_REL_GOOD = strprintf("../data/CYCLO_LLL_rel_good_%s_%s_new", deg_eq, size_roots);
  FILE_GP_GOOD  = strprintf("../data/CYCLO_GP_good_%s_%s_new", deg_eq, size_roots);

  FILE_LLL_ABS_BAD = strprintf("../data/CYCLO_LLL_abs_bad_%s_%s_new", deg_eq, size_roots);
  FILE_LLL_REL_BAD = strprintf("../data/CYCLO_LLL_rel_bad_%s_%s_new", deg_eq, size_roots);
  FILE_GP_BAD  = strprintf("../data/CYCLO_GP_bad_%s_%s_new", deg_eq, size_roots);

  print("deg_eq is: ", deg_eq);
  for (cond = conds[1], conds[2],		/* potential conductors */
	 print("***********************************");
       print("Conductor is ", cond);
       dim = eulerphi(cond);
       F = factor(cond);
       prime = F[1,1];
       dime = dim/eulerphi(cond/prime);
       print("Dimension is ", dimL);


       if ((deg_eq^(dime) > min(2^27, dime^4)), prime=1); /* avoid too large enumerations */
       
       if ((7 < dim) && (dim < dim_max),
	   
	   print(prime);
	   time_lll_abs = 0;
	   time_lll_rel = 0;
	   time_gp = 0;
	   
	   c_abs = 0;
	   c_rel = 0;

	   c_good = 0;
	   c_bad = 0;

	   pols = Cyclo_real_rel_ext_creation(cond, prime, 1, varg=z, vare=y, varK=t, varL=T);
	   print(pols);
	   pol_vec = pols[1..2];
	   p = pols[3];
	   

	   /* determine if Q(zeta_m) is a "good" field, i.e. it is cyclic */
	   my(b_good = cf_can_inert(cond););

	   for(i = 1, number_tests,
		 
		 /* ***************  first, relative setting  *************** */
		 if (deg_eq==2 || version=="split", /* split equation */
		     equation = Rel_equation_creation_split(pol_vec, deg_eq, size_roots)[1];,
		     
		     equation = Rel_equation_creation_split(pol_vec, 1, size_roots)[1]; /* equation w. only one root */
		     equation = equation*Rel_equation_creation_multiquad(pol_vec, (deg_eq-1)\2, size_roots)[1];
		     );
	     
	     
	       my(s=getabstime());
	       S_rel = Rel_pol_roots(pol_vec, equation, 0, 0);
	       time_lll_rel += (getabstime()-s)/1000.0;

	       
	       /* *************** then, absolute computation *************** */
	       p = subst(p, T, y);
	       
	       if (deg_eq==2 || version=="split",
		   equation = Equation_creation_nf_split(p, deg_eq, size_roots); ,

		   equation = Equation_creation_nf_split(p, 1, size_roots);
		   equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-1)\2, size_roots);
		   );
	     
	       /* first lll abs. */
	       my(s=getabstime());
	       S_abs = Solve_equation(p, equation, 0);
	       time_lll_abs += (getabstime()-s)/1000.0;
	     
	     
	       /* then gp (abs rep.) */
	       my(s=getabstime());
	       S = nfroots(p, equation);
	       time_gp += (getabstime()-s)/1000.0;
	     
	       if (vecprod(S)==vecprod(S_abs), c_abs += 1; );
	       if (length(S)==length(S_rel), c_rel += 1; );
	       kill(S);
	       kill(S_abs);
	       kill(S_rel);

	       
	       if (i%5==0, print(i "th element");
		   printf("Retrieved? %d \t %d \n", c_rel, c_abs);
		   printf("Times (gp / abs / rel) are: %s \t %s \t %s \n", time_gp/i, time_lll_abs/i, time_lll_rel/i);
		   );
	       );

	   time_gp /= number_tests;
	   time_lll_abs /= number_tests;
	   time_lll_rel /= number_tests;
	 
	   default(realprecision, 7);
	   
	   /* str_rel = strprintf("%s\t%s\t%s\t%s", cond, dim, time_lll_rel, c_rel); */
	   /* str_abs = strprintf("%s\t%s\t%s\t%s", cond, dim, time_lll_abs, c_abs); */
	   /* str_gp = strprintf("%s\t%s\t%s\t%s\t%s", cond, dim, time_gp, D_gp[1], D_gp[2]); */

	   str_rel = strprintf("%s\t%s\t%s\t%s", cond, dim, time_lll_rel, c_rel);
	   str_abs = strprintf("%s\t%s\t%s\t%s", cond, dim, time_lll_abs, c_abs);
	   str_gp = strprintf("%s\t%s\t%s", cond, dim, time_gp);

	   if (b_good,
	       f = fileopen(FILE_LLL_REL_GOOD, "a"); /* printing for heur method */
	       filewrite(f, str_rel);
	       fileclose(f);
	       
	       f = fileopen(FILE_LLL_ABS_GOOD, "a"); /* printing for heur method */
	       filewrite(f, str_abs);
	       fileclose(f);
	       
	       f = fileopen(FILE_GP_GOOD, "a"); /* printing for heur method */
	       filewrite(f, str_gp);
	       fileclose(f);  ,

	       f = fileopen(FILE_LLL_REL_BAD, "a"); /* printing for heur method */
	       filewrite(f, str_rel);
	       fileclose(f);
	       
	       f = fileopen(FILE_LLL_ABS_BAD, "a"); /* printing for heur method */
	       filewrite(f, str_abs);
	       fileclose(f);
	       
	       f = fileopen(FILE_GP_BAD, "a"); /* printing for heur method */
	       filewrite(f, str_gp);
	       fileclose(f);
	       
	       );

	   );
       );
};


/* L = Q(zeta_m) with m = prime^e for e in N  */
ExperimentsPrimePowerCyclotomics(prime, deg_eq, size_roots, number_tests, \
				 {dim_max=500}) = {
  my(p, pols, pol_vec, S, S_abs, S_rel, equation, field_dim, FILE_LLL_ABS, \
     FILE_LLL_REL, FILE_GP, time_lll_abs, time_lll_rel, time_gp, c_rel, c_abs, \
     times_gp, times_lll_rel, times_lll_abs);
  
  FILE_LLL_ABS_GOOD = strprintf("../data/PRIMEPOWER_CYCLO_LLL_abs_good_%s_%s_%s", \
				prime, deg_eq, size_roots);
  FILE_GP_GOOD  = strprintf("../data/PRIMEPOWER_CYCLO_GP_good_%s_%s_%s", prime,\
			    deg_eq, size_roots); 

  FILE_LLL_ABS_BAD = strprintf("../data/PRIMEPOWER_CYCLO_LLL_abs_bad_%s_%s_%s",\
			       prime, deg_eq, size_roots);
  FILE_GP_BAD  = strprintf("../data/PRIMEPOWER_CYCLO_GP_bad_%s_%s_%s", prime, \
			   deg_eq, size_roots); 

  exp_g_m = floor(log(10)/log(prime));
  
  exp_g = exp_g_m;
  exp_a = exp_g_m+1;
  dg = eulerphi(prime^exp_g_m);
  exp_e = exp_a - exp_g_m;
  de = prime^exp_e;
  da = dg*de;
  
  while (da < dim_max+1,
	 
	 FILE_LLL_REL_GOOD = strprintf("../data/PRIMEPOWER_CYCLO_LLL_rel_good_%s_%s_%s_%s",  prime, exp_e, deg_eq, size_roots);

	 FILE_LLL_REL_BAD = strprintf("../data/PRIMEPOWER_CYCLO_LLL_rel_bad_%s_%s_%s_%s", prime, exp_e, deg_eq, size_roots);

	 print("Dimension is ", da);
	 
	 time_lll_abs = 0;
	 time_lll_rel = 0;
	 time_gp = 0;
	 
 	 c_abs = 0;
	 c_rel = 0;
	 
	 pols = Cyclo_rel_ext_creation(prime, exp_g, exp_e);
	 print(pols);
	 pol_vec = pols[1..2];
	 p = pols[3];

	 /* determine if p is a "good" field for this equation size */
	 my(b_good = is_good_field(p, deg_eq)[1]);

	 
	 for(i = 1, number_tests,
	       
	       if (i%2==0,
		   printf ("%d th test\n", i-1);
		   printf("Time for rel. lll: %s\n", time_lll_rel/(i-1));
		   printf("Time for abs. lll: %s\n", time_lll_abs/(i-1));
		   printf("Time for gp: %s\n", time_gp/(i-1));
		   );
	     
	     if (deg_eq==2 || deg_eq!=prime,
		 equation = Rel_equation_creation_split(pol_vec, deg_eq, size_roots)[1];,
		 equation = Rel_equation_creation_split(pol_vec, 1, size_roots)[1];
		 equation = equation*Rel_equation_creation_multiquad(pol_vec, (deg_eq-1)\2, size_roots)[1];
		 );	   
	     
	     
	     my(s=getabstime());
	     S_rel = Rel_pol_roots(pol_vec, equation, 0, 0);
	     time_lll_rel += (getabstime()-s)/1000.0;
	     
	     
	     /* then absolute computation */
	     p = subst(p, x, y);

	     if (deg_eq==2 || prime!=deg_eq,
		 equation = Equation_creation_nf_split(p, deg_eq, size_roots); ,
		 
		 equation = Equation_creation_nf_split(p, 1, size_roots);
		 equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-1)\2, size_roots);
		 );
	     
	     /* first lll abs. */
	     my(s=getabstime());
	     S_abs = Solve_equation(p, equation, 0);
	     time_lll_abs += (getabstime()-s)/1000.0;
	     
	     
	     /* then gp (abs rep.) */
	     my(s=getabstime());
	     S = nfroots(p, equation);
	     time_gp += (getabstime()-s)/1000.0;
	     
	     if (vecprod(S)==vecprod(S_abs), c_abs += 1; );
	     if (length(S)==length(S_rel), c_rel += 1; );
	     kill(S);
	     kill(S_abs);
	     kill(S_rel);

	     );
	 
	 time_lll_rel /=  number_tests;
	 c_rel /=  1.0*number_tests;
	 
	 time_lll_abs /=  number_tests;
	 c_abs /=  1.0*number_tests;

	 time_gp /=  number_tests;
	 
	 default(realprecision, 7);
	 
	 /* str_rel = strprintf("%s\t%s\t%s\t%s", dg, de, time_lll_rel, c_rel); */
	 /* str_abs = strprintf("%s\t%s\t%s\t%s", dg, de, time_lll_abs, c_abs); */
	 /* str_gp = strprintf("%s\t%s\t%s\t%s\t%s", dg, de, time_gp, D_gp[1], D_gp[2]); */

	 str_rel = strprintf("%s\t%s\t%s\t%s", dg, de, time_lll_rel, c_rel);
	 str_abs = strprintf("%s\t%s\t%s\t%s", dg, de, time_lll_abs, c_abs);
	 str_gp = strprintf("%s\t%s\t%s", dg, de, time_gp);

	 
	 /* f = fileopen(FILE_LLL_REL, "a"); /\* printing for rel method *\/ */
	 /* filewrite(f, str_rel); */
	 /* fileclose(f); */
	 
	 /* f = fileopen(FILE_LLL_ABS, "a"); /\* printing for abs method *\/ */
	 /* filewrite(f, str_abs); */
	 /* fileclose(f); */

	 /* f = fileopen(FILE_GP, "a"); /\* printing for gp method *\/ */
	 /* filewrite(f, str_gp); */
	 /* fileclose(f); */


	 if (b_good, 
	     f = fileopen(FILE_LLL_REL_GOOD, "a"); /* printing for heur method */
	     filewrite(f, str_rel);
	     fileclose(f);
	       
	     f = fileopen(FILE_LLL_ABS_GOOD, "a"); /* printing for heur method */
	     filewrite(f, str_abs);
	     fileclose(f);
	       
	     f = fileopen(FILE_GP_GOOD, "a"); /* printing for heur method */
	     filewrite(f, str_gp);
	     fileclose(f);  ,

	     f = fileopen(FILE_LLL_REL_BAD, "a"); /* printing for heur method */
	     filewrite(f, str_rel);
	     fileclose(f);
	       
	     f = fileopen(FILE_LLL_ABS_BAD, "a"); /* printing for heur method */
	     filewrite(f, str_abs);
	     fileclose(f);
	       
	     f = fileopen(FILE_GP_BAD, "a"); /* printing for heur method */
	     filewrite(f, str_gp);
	     fileclose(f);  
	       
	     );

	 
	 while ((deg_eq^(de*prime) < 2^27) && (exp_g > exp_g_m),
		print("change ground field for relative method");
		/* update of the ground field */
		exp_g -= 1;
		dg = eulerphi(prime^exp_g);
		de = da/dg;
		exp_e = exp_a-exp_g;

		FILE_LLL_REL_GOOD = strprintf("../data/PRIMEPOWER_CYCLO_LLL_rel_good_%s_%s_%s_%s", prime, exp_e, deg_eq, size_roots);
		FILE_LLL_REL_BAD = strprintf("../data/PRIMEPOWER_CYCLO_LLL_rel_bad_%s_%s_%s_%s", prime, exp_e, deg_eq, size_roots);

		time_lll_rel = 0;
	 
		times_lll_rel = [];
	 
		c_rel = 0;
	 
		pols = Cyclo_rel_ext_creation(prime, exp_g, exp_e);
		pol_vec = pols[1..2];
		p = pols[3];
	 
		for(i = 1, number_tests,
	       
		      if (i%2==0,
			  printf ("%d th test\n", i-1);
			  printf("Time for rel. lll: %s\n", time_lll_rel/(i-1));
			  );

		    if (deg_eq==2 || deg_eq!=prime,
		    	equation = Rel_equation_creation_split(pol_vec, deg_eq, size_roots)[1];,
			equation = Rel_equation_creation_split(pol_vec, 1, size_roots)[1];
			equation = equation*Rel_equation_creation_multiquad(pol_vec, (deg_eq-1)\2, size_roots)[1];
			);

		    /* equation = Rel_equation_creation_split(pol_vec, deg_eq, size_roots)[1]; */
		    
		    my(s=getabstime());
		    S_rel = Rel_pol_roots(pol_vec, equation, 0, 0);
		    time_lll_rel += (getabstime()-s)/1000.0;
	     
		    if (length(S)==length(S_rel), c_rel += 1; );
		    kill(S_rel);
		    
		    );
		
		time_lll_rel /=  number_tests;
		c_rel /=  1.0*number_tests;	 
		
		default(realprecision, 7);
	 
		str_rel = strprintf("%s\t%s\t%s\t%s", dg, de, time_lll_rel, c_rel);
		
		if (b_good, 
		    f = fileopen(FILE_LLL_REL_GOOD, "a"); ,
		    f = fileopen(FILE_LLL_REL_BAD, "a");
		    );
		
		filewrite(f, str_rel);
		fileclose(f);
	
		);
	 exp_a += 1;
	 exp_e = 1;
	 exp_g = exp_a -1;
	 dg = eulerphi(prime^exp_g);
	 de = prime^exp_e;
	 da = dg*de;
	 );
  
};




/* ************************************************************************** */
/* ************************** KUMMER FIELDS ********************************  */
/* ************************************************************************** */

ExperimentsKummer_LLL_abs(exp, length_seq, deg_eq, vector_field, \
			  {number_tests=25, version="split"}) = {
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
	   
	   p = vector_field[i];
	   p = subst(p, z, y);

	   if (deg_eq==2 || nr==deg_eq,
	       equation = Equation_creation_nf_split(p, deg_eq, size_roots); ,

	       equation = Equation_creation_nf_split(p, nr, size_roots);
	       equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-nr)\2, size_roots);
	       /* equation = Equation_creation_nf_split(p, deg_eq, size_roots); , */
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



ExperimentsKummer_LLL_rel(exp, length_seq, deg_eq, vector_field, number_tests,\
			  {version="split"}) = {
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
	   pol_vec = vector_field[i];
	   pol_vec[2] = Mod(pol_vec[2], pol_vec[1]);
	   
	   if (deg_eq==2 || nr==deg_eq,
	       equation = Rel_equation_creation_split(pol_vec, deg_eq, size_roots)[1]; ,
	       equation = Rel_equation_creation_split(pol_vec, nr, size_roots)[1];
	       equation = equation*Rel_equation_creation_multiquad(pol_vec, (deg_eq-nr)\2, size_roots)[1];
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



ExperimentsKummer_GP(exp, length_seq, deg_eq, vector_field, number_tests,\
		     {version="split"}) = {
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
	   
	   p = vector_field[i];
	   p = subst(p, z, y);
	   
	   /* print(p); */
	   /* print(size_roots); */
	   if (deg_eq==2 || nr==deg_eq,
	       equation = Equation_creation_nf_split(p, deg_eq, size_roots); ,
	       
	       equation = Equation_creation_nf_split(p, nr, size_roots);
	       equation = equation*Equation_creation_nf_multiquad(p, (deg_eq-nr)\2, size_roots);
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
