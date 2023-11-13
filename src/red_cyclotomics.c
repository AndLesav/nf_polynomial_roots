/* ########################################################################## */
/* ############################### CYCLOTOMICS ############################## */
/* ########################################################################## */

/* simple versions of the functions for cyclotomic fields, 
   using the maximal real subfield K+ */

 /* creation of reduced basis for given precision */
cyclo_reduced_basis_latt(p, Binf, prec) =  {
  my(field_dim, B, M, C);
  field_dim = poldegree(p)\2;
  C = (field_dim-1)*(field_dim+2)\4;
  B = vector(field_dim, i, -round(Binf[i]<<prec));
  M = matconcat([B; C*matid(field_dim)]);
  M = qflll(M, 3);
  return(M);
};


/* Solve equation in QQ(zeta_m) */
CyclotomicsRelativeRoots(m, p, equation) = {
  my(prec, ns, v, field_dim, N, Ninv, zinf, einf, prec_add, L, Me, S);
  my(t = getabstime());
  /* p = polcyclo(m); */
  field_dim = poldegree(p);
  zeta_m = Mod(y, p);
  eta_m = zeta_m + 1/zeta_m;
  
  /* put eq in FF */
  [M, P] = FF_basis(p, 32, field_dim\4);
  pol1 = FF_pol_embedding_fam_nf(equation, M, P); 


  /* norm of transformation from basis (zeta_k + zeta^(-k))_k to (zeta^k)_k */
  N = matrix(field_dim, field_dim\2);
  for (i = 0, field_dim\2-1,
	 N[,i+1] = Vecrev(liftall(zeta_m^i + zeta_m^(-i) ), field_dim)~;
       );
  Ninv = N^(-1);

  /* print("inverse matrix computed in ", strtime(getabstime()-t)); */
  
  /* as always ; prec eval */
  [prec, ns, v] = round(Precision_eval(p, equation, 0));
  prec += round(log(log(field_dim\2))*log(field_dim\2)*(field_dim\2) + log(field_dim)/2*(field_dim\4));
  prec += round(log(norml2(Ninv))*field_dim);
  prec_add = v\2;
  
  default(realprecision, prec\2  + field_dim*field_dim  + prec_add);
  zinf = rootsof1(m)[2];
  einf = real(zinf+1/zinf);

  /* Bpinf = vector(field_dim\2, i, einf^(i-1)); */
  
  Bpinf = vector(field_dim\2, i, zinf^(i-1) + zinf^(1-i));
  Bpinf[1] = 1;
  Binf = vector(field_dim, i, zinf^(i-1));
  Bp = vector(field_dim\2, j, (zeta_m^(j-1) + zeta_m^(1-j)));

  Me = Mat([1, zinf; 1, conj(zinf)]);
  Meinv = Me^(-1);
  
  /* print("before reduction in ", strtime(getabstime()-t)); */
  L = cyclo_reduced_basis_latt(p, Bpinf, prec);
  /* print("reduction in ", strtime(getabstime()-t)); */
  
  QR = MyGSO_L2(L, 1, prec\2+field_dim);
   
  default(realprecision, round(prec\2+prec_add+field_dim^2));
  pol = Pol_embedding_nf(equation, Binf); /* put equation in RR */
  
  default(realprecision, round(prec\2)+ field_dim+prec_add\10);
  R = polroots(pol);        	    /* find roots in CC */
  /* print("before search in ", strtime(getabstime()-t)); */

    
  S = [];
  for(i=1, length(R),
	/* print("*******************************"); */
	
	my(mink = [R[i], conj(R[i])]); /* pick a root and its conj */

      ginf = Meinv*(mink~);


      default(realprecision, prec\3);
      /* deconding in K+ */
      [b1, g1] = MyDecodeBabai_2([real(ginf[1])], L, QR, prec, 1); /* first coord */
      g1 = g1~;

      if (b1,
	  /* second coord */
	  [b2, g2] = MyDecodeBabai_2([real(ginf[2])], L, QR, prec, 1); 
	  g2 = g2~;
	  );
      
      if (b1*b2==1,

	  g1 = g1[1] + sum(j=2, field_dim\2, g1[j]*Bp[j]);
	  g2 = g2[1] + sum(j=2, field_dim\2, g2[j]*Bp[j]);
	  
	  g = g1 + zeta_m*g2;
	  g1 = Vecrev(liftall(g), field_dim);

	  
	  gemb = FF_eval_pol(g1, pol1, M, P);
	  if(FF_test_vec(gemb, P),
	     S = (concat(S, [g]));
	     );
	  );
      if (number_solutions != 0,
	  if (length(S)==number_solutions, return(S););
	  );
      );      
  
  return(S);  
};
