/* FUNCTIONS TO SOLVE POLYNOMIAL EQUATIONS IN NUMBER FIELDS */

/* ABSOLUTE CASE */


/* Tries so solve an equation through an embedding in RR 
   prec = precision used
   basis = basis of nf(pol_field) in RR ( prec+dim^2 ) 
   lattice = [B | C*Id]~ with precision prec
   g = solution of eq */

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
Test_solve_equation_complex(pol_field, equation, prec, basis, lattice,  prec_add,\
			    {number_solutions}) = {
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


  





/* solving an equation in a number field with an embedding in RR 
   pol_field : polynomial defining the number field 
   equation : equation being solved
   version : 1 for certified or 0 for heuristic 
   number_solutions = number of solutions, if known */

Solve_equation_real(pol_field, equation, {version = 1, number_solutions=0}) = {
    local(prec, prec_add, n, latt, uni, base, field_dim, S, M, eq, S_full, ns, V, v);
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
    S = Test_solve_equation_real(pol_field, eq, prec, base, latt, prec_add, number_solutions, version);
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
    local(prec, prec_add, n, n0, latt, uni, base, base1, field_dim, S, S_full, ns, V, eq, v, n_t);
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
    
    /* print("precision:", prec); */
    prec_add = prec + (v \ 2);
    my(t = getabstime());
    [latt, uni, base, prec] = New_basis_complex(pol_field, prec, prec_add);
    /* print("Lattice created in: " strtime(getabstime()-t)); */
    
    S = Test_solve_equation_complex(pol_field, equation, prec, base, latt, prec_add, number_solutions);
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
}


/* ******************************************************************************** */
/* ******************************** RELATIVE  CASE ******************************** */
/* ******************************************************************************** */


Rel_pol_roots(pol_vec, eq, {method=1, cert_prec=1, ns=0}) = {
  my(dg, de, precision, B, M, r1, r1_g, lattice, uni, S, S_full, eq_n);
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
   
  my(s = getabstime());
  [precision, prec_add] = Rel_precision_eval(pol_vec, eq, cert_prec);
  if (!cert_prec, precision = precision + round(log(log(dg))*log(dg)*dg + dg*log(dg)/2)\1.5);
  /* print("precision computed in ", strtime(getabstime()-s)); */
  prec_add = prec_add;
  default(realprecision, round(precision));
  [B, M, r1, r1_g] = Rel_complex_minkowski(pol_vec, precision + prec_add + de*dg);
  /* print("mink computed in: ", strtime(getabstime()-s)); */

  uni = matid(dg);
  lattice = Rel_complex_basis_lattice(pol_vec, B, precision, uni, 1);
  /* print("lattice computed in: ", strtime(getabstime()-s)); */

  S = [];
  eq_n = eq;
  if(method==1,
     S_full = Rel_test_solve_equation_complex(pol_vec, eq, precision, prec_add, B, M, lattice[1], r1, ns);  ,
     S_full = Rel_test_solve_equation_complex_heur(pol_vec, eq, precision, prec_add, B, M, lattice[1], r1, r1_g, ns);
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
  R = Rel_pol_mink_roots(equation, pol_vec, Bg, Me, r1, precision\3+dg+
			 prec_add\10); /* ad hoc precision ? */

  dv = vector(length(R), i, length(R[i]));
  Me_inv = Me^(-1);
  default(realprecision, round(precision/3));
  S = [];
  old_e = vector(de, j, []);
   
  for(i=1, vecprod(dv),
	if(length(S)==ns && ns!=0, return(S));

      /* if(i%100==0, print(i " over " dp^de);); */
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
	  /* g1 = Rel_partial_decode(r, Me_inv, lattice, precision); */
	  [bt, g] = Rel_partial_decode_qr(r, Me_inv, lattice, QR, precision, 0);
	  /* g = Vec_to_multipol(g1, 2); */

	  gemb = FF_rel_eval_pol_fam(g, pol_vec, pol_emb, M, P);
	      
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
	t_a, QR, COUNT, dv, P, M, pol_emb);
    
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
  dp = poldegree(equation);
    
  my(s = getabstime());
  R = Rel_pol_mink_roots(equation, pol_vec, Bg, Me, r1, precision\3+dg+
			 prec_add\5); /* ad hoc precision ? */
  /* print("roots computation: ",  strtime(getabstime()-s)); */
    
  dv = vector(length(R), i, length(R[i]));
  Me_inv = Me^(-1);
    
  /* my(s = getabstime()); */
  QR = MyGSO_L2(lattice, 2, precision+dg);
  /* print("QR computation: ",  strtime(getabstime()-s)); */
    
  /* embedding in finite fields for verification */
  [M, P] = FF_rel_basis(pol_vec, 32, dg\4+1);
  pol_emb = FF_rel_pol_embedding_fam_nf(equation, pol_vec, M, P); /* put eq in FF */
    
  /* default(realprecision, round(precision\3)); */
  S = [];
  old_e = vector(de, j, [0]);
    
  for (i=1, vecprod(dv),
	 if(length(S)==ns && ns!=0, return(S));
	 
       /* if(i%10==0, print(i " over " vecprod(dv));); */
       /* e = IndexToExponent(i, dp, de); */
       e = IndexToExponent_multi(i, dv);
	 
       if (r1_g, e = NextConj_multi_simple(e, dv, r1); );
	 
       /* e = NextInd(e, dp, de, r1, old_e); */
       e = NextInd_multi(e, dv, old_e);
       
       /* i = sum(j=1, length(e), (e[j]-1)*dp^(de-j))+1; */
       i = ExponentToIndex_multi(e, dv);
       if (i > vecprod(dv), return(S));

       b=1;
       if (r1_g, b = TestExpoConj(e, dp, de, r1, old_e); );
       
       if (b,
	   
	   r = vector(de, j, R[j][e[j]]);
	   /* default(realprecision, round(precision\3)); */
	   
	   /* decode one coord. to check if this is a solution */
	   /* [bt, at, g] = Rel_test_one_coord(r, Me_inv, lattice, precision, nt, it, t_a); */
	   [bt, at, g] = Rel_test_one_coord_qr(r, Me_inv, lattice, QR, precision,  1);
	   

	   /* decode everything if a potential solution is detected */
	   if(at==1,
	      /* g1 = Rel_partial_decode(r, Me_inv, lattice, precision/\* , it, g *\/); */
	      [bt, g1] = Rel_partial_decode_qr(r, Me_inv, lattice, QR, precision, 1);

	      if (bt, 
		  g = Vec_to_multipol(g1, 2);
	      
		  /* verify if we got a true solution */
		  gemb = FF_rel_eval_pol_fam(g1, pol_vec, pol_emb, M, P);
		   
		  if(FF_test_vec(gemb, P),
		     S = (concat(S, [g]));
		     /* go to next potential index */
		     e_new = vector(length(e), j, 1);
		     e_new[1] = e[1]+1;
		     /* i = sum(j=1, length(e), (e_new[j]-1)*dp^(de-j)); */
		     i = ExponentToIndex_multi(e_new, dv)-1;
		      
		     for(j=1, de,
			   old_e[j] = vecsort(matconcat([old_e[j], e[j]])[1,]);
			 );,
		      
		     i = ExponentToIndex_multi(e, dv);
		     /* i = sum(j=1, length(e), (e[j]-1)*dp^(de-j))+1; */
		 
		     );
		  );
	      );
	   );
       );
  return(S);
};
