/* FUNCTIONS TO SOLVE POLYNOMIAL EQUATIONS IN NUMBER FIELDS */

/* ABSOLUTE CASE */


/* Tries so solve an equation through an embedding in RR 
   prec = precision used
   basis = basis of nf(pol_field) in RR ( prec+dim^2 ) 
   lattice = [B | C*Id]~ with precision prec
   g = solution of eq
*/

Test_solve_equation_real(pol_field, equation, prec, basis, lattice, prec_add, \
			 {number_solutions=0, version})=  {
  local(q, r, g1, gemb, pol, field_dim, R, S, M, P, pol1, m, G2);
  field_dim = poldegree(pol_field);
    
  my(s = getabstime());
  G2 = MyGSO_L2(lattice, 1, prec\2+2*field_dim);
  /* print("Rec GSO computation: ",  strtime(getabstime()-s)); */
    
  default(realprecision, round(prec\2+prec_add+field_dim^2));
  my(s = getabstime());
  pol = Pol_embedding_nf(equation, basis); /* put equation in RR */
  /* print("embedding computed in: ", strtime(getabstime()-s)); */

  [M, P] = FF_basis(pol_field, 32, field_dim\4);
  pol1 = FF_pol_embedding_fam_nf(equation, M, P); /* put eq in FF */
         
  default(realprecision, round(prec\3)+ field_dim+prec_add\10);
  R = polrootsreal(pol);        	    /* find roots in RR */
  /* print("embedding and roots computed in: ", strtime(getabstime()-s)); */
  S= [];
  
  for(i=1, length(R),
	my(s = getabstime());
      [b, g1] = MyDecodeBabai_2([R[i]], lattice, G2, prec, !version);
      g1 = g1~;
      /* print("Decoding with GSO done in: ", strtime(getabstime()-s)); */

      gemb = FF_eval_pol(g1, pol1, M, P);
      g1 = Mod(Polrev(g1, y), pol_field);
              
      my(s = getabstime());
      if(FF_test_vec(gemb, P),
	 /* print("verification done in: ", strtime(getabstime()-s)); */
	 S = (concat(S, [g1]));
	 );
      if (number_solutions != 0,
	  if (length(S)==number_solutions, return(S););
	  );
      );
  return(S);
};


/* save as before but embedding in CC */
Test_solve_equation_complex(pol_field, equation, prec, basis, lattice,\
			    prec_add, {number_solutions}) = {
  local(g1, pol, field_dim, R, S, gemb, pol1, M, P, b, G2);
  field_dim = poldegree(pol_field);

  [M, P] = FF_basis(pol_field, 32, field_dim\4);
  pol1 = FF_pol_embedding_fam_nf(equation, M, P); /* put eq in FF */
        
  my(s = getabstime());
  G2 = MyGSO_L2(lattice, 2, prec+field_dim);
  /* print("Rec GSO computation: ",  strtime(getabstime()-s)); */

  default(realprecision, round(prec\3+prec_add+field_dim^2));
  my(s = getabstime());
  pol = Pol_embedding_nf(equation, basis); /* put equation in CC */
  /* print("embedding computed in: ", strtime(getabstime()-s)); */
	
  default(realprecision, round(prec\3)+field_dim + prec_add\10);
  R = polroots(pol);		    /* find roots in CC */
  /* print("embedding and roots computed in: ", strtime(getabstime()-s)); */
  
  S= [];
  b = 1;
  for(i=1, length(R),
	my(s = getabstime());
      [b, g1] = MyDecodeBabai_2([real(R[i]), imag(R[i])], lattice, G2, prec);
      g1 = g1~;
      /* print("Decoding with GSO done in: ", strtime(getabstime()-s)); */
	
      if (b, 
	  my(s = getabstime());
	  gemb = FF_eval_pol(g1, pol1, M, P);
	  g1 = Mod(Polrev(g1, y), pol_field);
	  if(FF_test_vec(gemb, P),
	     /* print("verification done in: ", strtime(getabstime()-s)); */
	     S = (concat(S, [g1]));
	     );
	  if (number_solutions != 0,
	      if (length(S)==number_solutions, return(S););
	      );
	  );
      );
  return(S);
};


/* ##########################################################################  */
/*                     FUNCTIONS FOR TESTS REGARDING EARLY ABORT               */
/* ##########################################################################  */


/*
  FUNCTION FOR EXPERIMENTS RELATED TO EARLY ABORT
  Tries so solve an equation through an embedding in RR 
  prec = precision used
  basis = basis of nf(pol_field) in RR ( prec+dim^2 ) 
  lattice = [B | C*Id]~ with precision prec
  g = solution of eq
  return as well the max number of steps in nearest plane when wpplied to
  complex roots not corresponding to any element :
*/

Test_solve_equation_real_nbsteps(pol_field, equation, prec, basis, lattice,\
				 prec_add, {number_solutions=0, version})=  {
  local(q, r, g1, gemb, pol, field_dim, R, S, M, P, pol1, m, G2);
  field_dim = poldegree(pol_field);
    
  my(s = getabstime());
  G2 = MyGSO_L2(lattice, 1, prec\2+2*field_dim);
  /* print("Rec GSO computation: ",  strtime(getabstime()-s)); */
    
  default(realprecision, round(prec\2+prec_add+field_dim^2));
  my(s = getabstime());
  pol = Pol_embedding_nf(equation, basis); /* put equation in RR */
  /* print("embedding computed in: ", strtime(getabstime()-s)); */

  [M, P] = FF_basis(pol_field, 32, field_dim\4);
  pol1 = FF_pol_embedding_fam_nf(equation, M, P); /* put eq in FF */
         
  default(realprecision, round(prec\3)+ field_dim+prec_add\10);
  R = polrootsreal(pol);        	    /* find roots in RR */
  /* print("embedding and roots computed in: ", strtime(getabstime()-s)); */
  S= [];
  nb_steps = 0;
  for(i=1, length(R),
	my(s = getabstime());
      [b, g1, nb_s] = MyDecodeBabai_nbsteps([R[i]], lattice, G2, prec, !version);
      g1 = g1~;
      /* print("Decoding with GSO done in: ", strtime(getabstime()-s)); */

      gemb = FF_eval_pol(g1, pol1, M, P);
      g1 = Mod(Polrev(g1, y), pol_field);
              
      my(s = getabstime());
      if(FF_test_vec(gemb, P),
	 /* print("verification done in: ", strtime(getabstime()-s)); */
	 S = concat(S, [g1]);,
	 nb_steps = max(nb_steps, nb_s);
	 );
      if (number_solutions != 0,
	  if (length(S)==number_solutions, return(S););
	  );
      );
  return([S, nb_steps]);
};

/*  FUNCTION TEST EARLY ABORT
  save as before but embedding in CC */
Test_solve_equation_complex_nbsteps(pol_field, equation, prec, basis, lattice, \
				    prec_add, {number_solutions, version}) = {
  local(g1, pol, field_dim, R, S, gemb, pol1, M, P, b, G2);
  field_dim = poldegree(pol_field);

  [M, P] = FF_basis(pol_field, 32, field_dim\4);
  pol1 = FF_pol_embedding_fam_nf(equation, M, P); /* put eq in FF */
        
  my(s = getabstime());
  G2 = MyGSO_L2(lattice, 2, prec+field_dim);
  /* print("Rec GSO computation: ",  strtime(getabstime()-s)); */

  default(realprecision, round(prec\3+prec_add+field_dim^2));
  my(s = getabstime());
  pol = Pol_embedding_nf(equation, basis); /* put equation in CC */
  /* print("embedding computed in: ", strtime(getabstime()-s)); */
	
  default(realprecision, round(prec\3)+field_dim + prec_add\10);
  R = polroots(pol);		    /* find roots in CC */
  /* print("embedding and roots computed in: ", strtime(getabstime()-s)); */
  
  S=[];
  b = 1;
  nb_steps = 0;
  for(i=1, length(R),
	my(s = getabstime());
      [b, g1, nb_s] = MyDecodeBabai_nbsteps([real(R[i]), imag(R[i])], lattice,\
					    G2, prec, !version);
      g1 = g1~;
      /* print("Decoding with GSO done in: ", strtime(getabstime()-s)); */
	
      my(s = getabstime());
      gemb = FF_eval_pol(g1, pol1, M, P);
	  g1 = Mod(Polrev(g1, y), pol_field);
      if(FF_test_vec(gemb, P),
	 /* print("verification done in: ", strtime(getabstime()-s)); */
	 S = (concat(S, [g1]));,
	 nb_steps = max(nb_steps, nb_s);
	 );
      if (number_solutions != 0,
	  if (length(S)==number_solutions, return(S););
	  );
      );
  return([S, nb_steps]);
};
  
/* ########################################################################  */
/* ########################################################################  */
  
  
  
  
  
  
/* solving an equation in a number field with an embedding in RR 
   pol_field : polynomial defining the number field 
   equation : equation being solved
   version : 1 for certified or 0 for heuristic 
   number_solutions = number of solutions, if known */

Solve_equation_real(pol_field, equation, {version = 1, number_solutions=0}) = {
  local(prec, prec_add, n, latt, uni, base, field_dim, S, M, eq, S_full, ns, V,
	v);
  my(t = getabstime());
  field_dim = poldegree(pol_field);
  n = 0;
  [prec, ns, v] = round(Precision_eval(pol_field, equation, version));
  /* print("precision evaluated in: ", , strtime(getabstime()-t)); */
  if (!version,
      prec = prec + ceil(log(log(field_dim))*log(field_dim)*field_dim
			 + (field_dim\2)*log(field_dim)\2
			 ); 
      );
  prec_add = prec + v\2;
        
  my(t = getabstime());   
  [latt, uni, base, prec] = New_basis(pol_field, prec, prec_add);
  /* print("Lattice created in: ", strtime(getabstime()-t)); */
  
  eq = equation;
  S = Test_solve_equation_real(pol_field, eq, prec, base, latt, prec_add,
			       number_solutions, version);
  /* print("equation maybe solved in: " strtime(getabstime()-t)); */
  S_full = S;
   
  if (number_solutions != 0,
      while (length(S_full) < number_solutions,
	     eq = Update_equation(pol_field, eq, S);
	     [latt, uni, base, prec] = Real_basis_lattice_update_geom(pol_field, prec, prec_add, sqrt(exp(1)), 1, uni, 0);
	     S = Test_solve_equation_real(pol_field, eq, prec, base, latt, prec_add, number_solutions-length(S_full));
	     S_full = concat(S_full, S);
	     ); 
      );
  /* print("roots computed in: " strtime(getabstime()-t) "\n"); */
  return(S_full);
};



/* same as before but with embedding in CC \ RR */
Solve_equation_complex(pol_field, equation, {version=1, number_solutions=0})=
  {
    local(prec, prec_add, n, n0, latt, uni, base, base1, field_dim, S, S_full,
	  ns, V, eq, v, n_t);
    field_dim = poldegree(pol_field);
    deg = poldegree(equation);
    eq = equation;
    
    my(t = getabstime());
    [prec, ns, v] = Precision_eval(pol_field, equation, version);
    /* print("precision computed in: " strtime(getabstime()-t)); */

    if (!version,
	prec = prec\2 + ceil(log(log(field_dim))*log(field_dim)*field_dim
			     + field_dim * log(field_dim)/2
			     );
	);
    
    prec_add = prec + (v \ 2);
    my(t = getabstime());
    [latt, uni, base, prec] = New_basis_complex(pol_field, prec, prec_add);
    /* print("Lattice created in: " strtime(getabstime()-t)); */
    
    S = Test_solve_equation_complex(pol_field, equation, prec, base, latt,\
				    prec_add, number_solutions);
    /* print("roots computed in: " strtime(getabstime()-t)); */

    S_full = S;
    if (number_solutions != 0,
	while (length(S_full) < number_solutions,
	       eq = Update_equation(pol_field, eq, S);
	       [latt, uni, base, prec] = Complex_basis_lattice_update_geom(pol_field, prec, prec_add, sqrt(2), 1, uni);
	       S = Test_solve_equation_complex(pol_field, eq, prec, base, latt, prec_add, number_solutions - length(S_full), base1);
	       S_full = concat(S_full, S);
	       );
	);
    /* print("roots computed in: " strtime(getabstime()-t)); */
    return(S_full);
  };




/* solve polynomial equation in number field 
   if r1>0 pick a real embedding
*/
Solve_equation(pol_field, equation, {version=1, number_solutions=0}) = {
  r = MySignature(pol_field);
  if (r[1]==0,
      return(Solve_equation_complex(pol_field, equation, version, number_solutions)); ,
      return(Solve_equation_real(pol_field, equation, version, number_solutions));,
      );
};


/* ************************************************************************** */
/* **************************** RELATIVE  CASE ****************************** */
/* ************************************************************************** */


/* solve polynomial equation eq in relative extension of number fields
   defined by pol_vec
   method : 0 (heuristic) or 1 (certified)
   cert_prec : 0 or 1 (certified)
   ns : if not 0, maximum number of solutions to find */
Rel_pol_roots(pol_vec, eq, {method=1, cert_prec=1, ns=0}) = {
  my(dg, de, precision, B, M, r1, r1_g, lattice, uni, S, S_full, eq_n);
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
  
  my(s = getabstime());
  [precision, prec_add] = Rel_precision_eval(pol_vec, eq, cert_prec);
  
  if (!cert_prec,
      precision = precision + round(log(log(dg))*log(dg)*dg/2 + dg*log(dg)/2);
      precision = round(1.15*precision);
      precision = max(precision, 200);
      );
  /* print("precision computed in ", strtime(getabstime()-s)); */
  /* print("precision is ", precision); */
  prec_add *= 2;
  prec_add = round(prec_add);

  default(realprecision, round(precision));
  [B, M, r1, r1_g] = Rel_complex_minkowski(pol_vec, precision + prec_add + 2*de*dg);
  /* print("mink computed in: ", strtime(getabstime()-s)); */

  uni = matid(dg);
  lattice = Rel_complex_basis_lattice(pol_vec, B, precision, uni, 1);
  /* printf("lattice computed in: %s \n\n", strtime(getabstime()-s)); */
  
  S = [];

  eq_n = eq;

  if(method==1,
     /* S_full = Rel_test_solve_equation_complex(pol_vec, eq, precision, prec_add, B, M, lattice[1], r1, ns);  , */
     S_full = Rel_test_solve_equation_complex_new(pol_vec, eq, precision, prec_add, B, M, lattice[1], r1, ns);  ,
     
     S_full = Rel_test_solve_equation_complex_heur_new(pol_vec, eq, precision, prec_add, B, M, lattice[1], r1, r1_g, ns);
     /* S_full = Rel_test_solve_equation_complex_heur(pol_vec, eq, precision, prec_add, B, M, lattice[1], r1, r1_g, ns); */
     );
  
  if (ns != 0,
      while (length(S_full)!=ns,
	     eq_n = Update_equation(pol_field, eq_n, S);
	     precision = round(1.25*precision);
	     default(realprecision, round(precision));
	     [B, M, r1] = Rel_complex_minkowski(pol_vec, round(5*precision));
	     lattice = Rel_complex_basis_lattice(pol_vec, B, precision, uni);
	     if(method==1,
		S = Rel_test_solve_equation_complex(pol_vec, eq, precision, B, M, lattice[1], r1, ns-length(S_full));  ,
		S = Rel_test_solve_equation_complex_heur(pol_vec, eq, precision, B, M, lattice[1], r1, ns-length(S_full));
		);
	     S_full = concat(S_full, S);
	     );
      );
  return(S_full);
};


/* Try to solve equation in relative extension of number fields L/K : certified
   pol_vec : [PK, PL] PK defining pol of K/QQ and PL defining pol of L/K
   equation : equation to solve
   precision : precision at which computations are done
   prec_add : additional prec used to ensure correctness
   Bg : basis of ground field K
   Me : relative minkowski
   lattice : basis lattice used to decode
   r1 : r1(L)
   r1_g : r1(K)
*/
Rel_test_solve_equation_complex(pol_vec, equation, precision, prec_add, Bg, Me,\
				lattice, r1, r1_g, {ns=poldegree(equation)}) = {
  local(dg, de, dp, R, S, e, r, Me_inv, g, b, old_e, M, P, g_emb, g1);
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
  dp = poldegree(equation);


  /* embedding in finite fields for verification */
  [M, P] = FF_rel_basis(pol_vec, 32, dg*de\4);
  pol_emb = FF_rel_pol_embedding_fam_nf(equation, pol_vec, M, P); /* put eq in FF */

  /* QR decomp for babai nearest plane */
  QR = MyGSO_L2(lattice, 2, precision\2+dg);
    
  /* default(realprecision, round(precision\3)); */
  R = Rel_pol_mink_roots(equation, pol_vec, Bg, Me, r1, precision\2+dg+
			 prec_add\5); /* ad hoc precision ? */

  dv = vector(length(R), i, length(R[i]));
  Me_inv = Me^(-1);
  default(realprecision, round(precision/3));
  S = [];
  old_e = vector(de, j, []);
   
  for(i=1, vecprod(dv),
	if(length(S)==ns && ns!=0, return(S));

      /* if(i%1000==0, print(i " over " dp^de);); */
      /* print(e); */
      e = IndexToExponent_multi(i, dv);
      /* e = IndexToExponent(i, dp, de); */
      if (r1_g, e = NextConj_multi_simple(e, dv, r1); );

      e = NextInd_multi(e, dv, old_e);
      /* e = NextInd(e, dp, de, r1, old_e); */

      i = ExponentToIndex_multi(e, dv);
      /* i = sum(j=1, length(e), (e[j]-1)*dp^(de-j))+1; */

      if (i > vecprod(dv), return(S));

      b=1;
      if(r1_g != 0, b = TestExpoConj(e, dp, de, r1, old_e););
      if (b,
	  r = vector(de, j, R[j][e[j]]);

	  default(realprecision, round(precision/3));

	  g1 = Rel_partial_decode(r, Me_inv, lattice, precision);
	  /* [bt, g] = Rel_partial_decode_qr(r, Me_inv, lattice, QR, precision, 0); */
	  /* g = Vec_to_multipol(g1, 2); */

	  gemb = FF_rel_eval_pol_fam(g1, pol_vec, pol_emb, M, P);
	  
	  if(FF_test_vec(gemb, P),
	     /* if(subst(equation, x, g)==0, */
	     /* input(); */
	     /* print("TRUE"); */
	     S = (concat(S, [g]));
	     e_new = vector(length(e), j, 1);
	     e_new[1] = e[1]+1;
	     i = ExponentToIndex_multi(e_new, dv)-1;
	     /* i = sum(j=1, length(e), (e_new[j]-1)*dp^(de-j)); */
	     for(j=1, de,
		   old_e[j] = matconcat([old_e[j], e[j]])[1,];
		 /* print(old_e); */
		 );,
	     /* i = sum(j=1, length(e), (e[j]-1)*dp^(de-j))+1; */
	     i = ExponentToIndex_multi(e, dv);
	     );
	   
	  /* print(length(S)); */
	  if(length(S)==ns && ns!=0, return(S));
	  );
      );
  return(S);
};



  

Rel_test_solve_equation_complex_heur(pol_vec, equation, precision, prec_add, Bg, Me,\
				     lattice, r1, r1_g, {ns=poldegree(equation)}) = {
  local(dg, de, dp, R, S, e, r, Me_inv, g, gemb, g1, b, old_e, nt, it, bt, at, i_new,\
	t_a, QR, COUNT, dv, P, M, pol_emb, Z_inv);
    
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
  dp = poldegree(equation);
  
  my(s = getabstime());

  R = Rel_pol_mink_roots(equation, pol_vec, Bg, Me, r1, precision\3+dg+
			 prec_add\10); /* ad hoc precision ? */
  /* print("roots computation: ",  strtime(getabstime()-s)); */


  dv = vector(length(R), i, length(R[i]));
  pv = myProdvec(dv);
  
  Me_inv = Me^(-1);
  
  Z_inv = BatchMatMult(Me_inv, R);
  
  
  /* my(s = getabstime()); */
  QR = MyGSO_L2(lattice, 2, precision\2+dg);
  /* print("QR computation: ",  strtime(getabstime()-s)); */

  /* embedding in finite fields for verification */
  [M, P] = FF_rel_basis(pol_vec, 32, dg\4+1);
  pol_emb = FF_rel_pol_embedding_fam_nf(equation, pol_vec, M, P); /* put eq in FF */

    
  /* default(realprecision, round(precision\3)); */
  S = [];
  old_e = vector(de, j, [0]);
  count = 0;

  my(myc=0 , myT=0);

  for (i=1, vecprod(dv),
	 
	 if(length(S)==ns && ns!=0, return(S));
       
       /* if(i%10^4==0, print(i " over " vecprod(dv));); */
   
       default(realprecision, 10);
       e = IndexToExponent_multi(i, dv, pv);
       if (r1_g, e = NextConj_multi_simple(e, dv, r1); );
       
       my(s = getabstime());
       e = NextInd_multi(e, dv, old_e);

       i = ExponentToIndex_multi(e, dv, pv);
       if (i > vecprod(dv), return(S));

       b=1;
       if (r1_g, b = TestExpoConj(e, dp, de, r1, old_e); );
       if (b,
	   
	   r = vector(de, j, R[j][e[j]]);
	   /* r = (vector(de, j, Z_inv[j][e[j]])); */
	   
	   /* decode one coord. to check if this is a solution */
	   [bt, at, g] = Rel_test_one_coord_qr(r, Me_inv, lattice, QR, precision,  1);
	   at = 1;
	   
	   /* decode everything if a potential solution is detected */
	   if(at==1,
	      my(s = getabstime());
	      [bt, g1] = Rel_partial_decode_qr(r, Me_inv, lattice, QR, precision, 1);
	      
	      if (bt,
		  
		  g = Vec_to_multipol(g1, 2);
		  count += 1;
		  /* printf("there is maybe a solution: %o\n", count); */
		  /* print("decoding done in: ",  strtime(getabstime()-s)); */
		  		  
		  my(s = getabstime());
		  /* verify if we got a true solution */
		  gemb = FF_rel_eval_pol_fam(g1, pol_vec, pol_emb, M, P);
		  
		  if(FF_test_vec(gemb, P),
		     
		     S = (concat(S, [g]));
		     /* go to next potential index */
		     e_new = vector(length(e), j, 1);
		     e_new[1] = e[1]+1;
		     /* i = sum(j=1, length(e), (e_new[j]-1)*dp^(de-j)); */
		     i = ExponentToIndex_multi(e_new, dv, pv)-1;
		     
		     for(j=1, de,
			   old_e[j] = vecsort(matconcat([old_e[j], e[j]])[1,]);
			 );,
		     
		     i = ExponentToIndex_multi(e, dv, pv);
		     );
		  /* print("verif done in: ",  strtime(getabstime()-s)); */
		  );
	      );
	   );
       );
  return(S);
};


Rel_first_sol(~R, lattice, QR, precision, M, P, pol_emb, dp, de, r1_g, r1) = {
  my(S, dv, pv);
  
  dv = vector(length(R), i, length(R[i]));
  /* dv = vector(length(Zinv), i, length(Zinv[i])); */

  pv = myProdvec(dv);
  pr = pv[#pv];
  e = vector(length(dv), i, 1);
  e[#e] = 0;

  for (i=1, pr,
	 
	 default(realprecision, 10);
       NextExp_multi(~e, dv);

       /* update */
       k = 2;
       bool = e[2]==1;
       while (bool && k < #e,
	      k += 1;
	      bool = bool && (e[k]==1);
	      );
	      
       if (bool && e[1] != 1,
	   listpop(~R[1], e[1]-1);
	   dv[1] -= 1;
	   if (r1_g,
	       listpop(~R[2], e[1]-1);
	       dv[2] -= 1;
	       );
	   pv = myProdvec(dv);
	   pr = pv[#pv];
	   e = vector(length(dv), i, 1);
	   );


       if (r1_g, e = NextConj_multi_simple(e, dv, r1);
	   i = ExponentToIndex_multi(e, dv, pv);
	   );
       
       /* if (i > vecprod(dv), return([0, S])); */

       b=1;
       
       /* if (r1_g, b = TestExpoConj(e, dp, de, r1); ); */

       if (b,
	   r = vector(de, j, R[j][e[j]]);
	   
	   my(s = getabstime());
	   [bt, g1] = Rel_partial_decode_qr(r, Me_inv, lattice, QR, precision, 0);
	   g = Vec_to_multipol(g1, 2);
	   /* print("decoding done in: ",  strtime(getabstime()-s)); */
	   
	   my(s = getabstime());
	   /* verify if we got a true solution */
	   gemb = FF_rel_eval_pol_fam(g1, pol_vec, pol_emb, M, P);
	   
	   if(FF_test_vec(gemb, P),
	      S = g;
	      for (j =1, length(R),
		     listpop(~R[j], e[j]);
		   );
	      return([1, S]);
	      );
	   );
       );
  return([0, S])
};


Rel_first_sol_heur(~R, lattice, QR, precision, M, P, pol_emb, dp, de, r1_g, r1) = {
  my(S, dv, pv, r, e, pr, bool, k);
  my(s_time = getabstime());
  count = 0;
  dv = vector(length(R), i, length(R[i]));
  pv = myProdvec(dv);
  pr = pv[#pv];
  e = vector(length(dv), i, 1);
  e[#e] = 0;
  for (i=1, pr,
	 
	 default(realprecision, 10);
       my(s = getabstime());
       /* e = IndexToExponent_multi(i, dv, pv); */
       NextExp_multi(~e, dv);

       /* update */
       k = 2;
       bool = e[2]==1;
       while (bool && k < #e,
	      k += 1;
	      bool = bool && (e[k]==1);
	      );
	      
       if (bool && e[1] != 1,
	   /* print(dv); */
	   listpop(~R[1], e[1]-1);
	   dv[1] -= 1;
	   if (r1_g && r1==0,
	       listpop(~R[2], e[1]-1);
	       dv[2] -= 1;
	       );
	   pv = myProdvec(dv);
	   pr = pv[#pv];
	   e = vector(length(dv), i, 1);
	   i = 1
	   );

       /* if (vecprod(e)==1 && i !=1,  return([0, S])); */
       if (i > pr, return([0, S]));

       if (r1_g,
	   e = NextConj_multi_simple(e, dv, r1);

	   i = ExponentToIndex_multi(e, dv, pv);
	   
	   );

       b=1;
       
       /* if (r1_g, b = TestExpoConj(e, dp, de, r1); ); */

       if (b,
	   r = vector(de, j, R[j][e[j]]);

	   s = getabstime();

	   /* decode one coord. to check if this is a solution */
	   [bt, at, g] = Rel_test_one_coord_qr(r, Me_inv, lattice, QR, precision, 1);
	   /* print("test one coord done in: ", strtime(getabstime() - s)); */
	   
	   /* decode everything if a potential solution is detected */
	   if(at==1,
	      count += 1;
	      /* printf("there is maybe a solution\n"); */
	      my(ss = getabstime());
	      my(s = getabstime());


	      [bt, g1] = Rel_partial_decode_qr(r, Me_inv, lattice, QR, precision, 1);
	      
	      if (bt,
		  g = Vec_to_multipol(g1, 2);
		  /* print("decoding done in: ",  strtime(getabstime()-s)); */
		  		  
		  my(s = getabstime());
		  /* verify if we got a true solution */
		  gemb = FF_rel_eval_pol_fam(g1, pol_vec, pol_emb, M, P);
		  /* print("embedding done in: ", strtime(getabstime()-s)); */

		  my(s = getabstime());

		  if(FF_test_vec(gemb, P),

		     /* print("one sol verified in FF: ", strtime(getabstime()-s)); */

		     /* print("one sol verified in : ", strtime(getabstime()-ss)); */

		     S = g;
		     for (j =1, length(R),
			    listpop(~R[j], e[j]);
			  
			  );

		     /* pv = myProdvec(dv); */
		     /* pr = pv[#pv]; */
		     /* e = vector(length(dv), i, 1); */
		     /* e[#e]=0; */
		     /* i = 0; */

		     /* print("end of work for one sol : ", strtime(getabstime()-ss)); */
		     return([1, S]);
		     );
		  );
	      );
	   );
       );
  return([0, S])
};


Rel_test_solve_equation_complex_new(pol_vec, equation, precision, prec_add, Bg,\
					 Me, lattice, r1, r1_g,\
					 {ns=poldegree(equation)}) = {
  
  local(dg, de, dp, R, S, e, r, Me_inv, g, gemb, g1, b, old_e, nt, it, bt, at,\
	i_new, t_a, QR, COUNT, dv, P, M, pol_emb, Z_inv);
    
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
  dp = poldegree(equation);

  /* embedding in finite fields for verification */
  [M, P] = FF_rel_basis(pol_vec, 32, dg\4+1);
  pol_emb = FF_rel_pol_embedding_fam_nf(equation, pol_vec, M, P); /* put eq in FF */

  /* my(s = getabstime()); */
  QR = MyGSO_L2(lattice, 2, precision\2+dg);

  my(s = getabstime());
  
  R = Rel_pol_mink_roots(equation, pol_vec, Bg, Me, r1, precision\3+dg+
			 prec_add\10); /* ad hoc precision ? */
  /* print("roots computation: ",  strtime(getabstime()-s)); */

  dv = vector(length(R), i, length(R[i]));
  pv = myProdvec(dv);

  Me_inv = Me^(-1);
  Z_inv = BatchMatMult(Me_inv, R);
  /* print("mink inv * roots computed in: ",  strtime(getabstime()-s)); */

  /* print("QR computation: ",  strtime(getabstime()-s)); */

  kill(Me_inv);
  kill(Me);
  kill(R);

    

  S = [];
  my(tb = 1);

  while ( (tb && ((length(S)<ns && ns!=0) || length(S)<dp)),
	  [tb, so] = Rel_first_sol(~Z_inv, lattice, QR, precision, M, P, pol_emb, \
				   dp, de, r1_g, r1);
	  S = concat(S, [so]);
	  );
   
  return(S);
};



Rel_test_solve_equation_complex_heur_new(pol_vec, equation, precision, prec_add, Bg,\
					 Me, lattice, r1, r1_g,\
					 {ns=poldegree(equation)}) = {
  local(dg, de, dp, R, S, e, r, Me_inv, g, gemb, g1, b, old_e, nt, it, bt, at, i_new,\
	t_a, QR, COUNT, dv, P, M, pol_emb, Z_inv);
    
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
  dp = poldegree(equation);
  
  /* embedding in finite fields for verification */
  [M, P] = FF_rel_basis(pol_vec, 32, 4);
  pol_emb = FF_rel_pol_embedding_fam_nf(equation, pol_vec, M, P); /* put eq in FF */

  /* my(s = getabstime()); */
  QR = MyGSO_L2(lattice, 2, precision\3+dg);
  /* print("QR computation: ",  strtime(getabstime()-s)); */
  
  my(s = getabstime());
  
  R = Rel_pol_mink_roots(equation, pol_vec, Bg, Me, r1, precision\2+dg+
			 prec_add\10); /* ad hoc precision ? */
  /* print("roots computation: ",  strtime(getabstime()-s)); */
  
  default(realprecision, 10);
  
  dv = vector(length(R), i, length(R[i]));
  pv = myProdvec(dv);
  Me_inv = Me^(-1);
  
  Z_inv = BatchMatMult(Me_inv, R);
  /* print("mink inv * roots computed in: ",  strtime(getabstime()-s)); */
    
  kill(Me_inv);
  kill(Me);
  kill(R);
  Z_inv = Vec_shuffle(Z_inv, r1_g, r1);
  /* print("shuffle computed in: ",  strtime(getabstime()-s)); */
    
  S = [];
  my(tb = 1);

  /* S = Rel_first_sol_heur(~Z_inv, lattice, QR, precision, M, P, pol_emb, \ */
  /* 			 dp, de, r1_g, r1); */

  while ((tb && ((length(S)<ns && ns!=0) || length(S)<dp)),
	  /* [tb, so] = Rel_first_sol_heur(~R, lattice, QR, precision, M, P, pol_emb, \ */
				   /* dp, de, r1_g, r1); */
	  [tb, so] = Rel_first_sol_heur(~Z_inv, lattice, QR, precision, M, P, pol_emb,\
					dp, de, r1_g, r1);
	  
	  if (tb,  S = concat(S, [so]));
	  
	  /* print(  vector(length(Z_inv), i, length(Z_inv[i]))); */
	  );
  return(S);
};
