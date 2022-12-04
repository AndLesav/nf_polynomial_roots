/* pol roots bound real */
MyRealPol_roots_bound_vec(pol_mink, r1)={
  my(v);
  v = vector(r1);
  for(i=1, r1,
	v[i] = vecmax(polrootsreal(pol_mink[i]));
      );
  return(v);
};


/* evaluation size root */
Size_root_eval_split(equation)=  {
  return(log(norml2(liftall(polcoeff(equation, poldegree(equation)-1))))/(log(2)));
};


/* bound on roots of vector of polynomials */
Polroots_bound_vec(pol_vec) = {
  /* default(realprecision, prec); */
  return(vector(length(pol_vec), i, polrootsbound(pol_vec[i])));
};


/* bound on roots of vector of polynomials */
Polroots_norm_vec(pol) = {
  my(R = polroots(pol));
  return(vecsort(vector(length(R), i, norml2(R[i]))));
};


/* bound on roots of vector of polynomials */
Polroots_norm_mat(pol_vec) = {
  my(M = Polroots_norm_vec(pol_vec[1]));
  for (i=2, length(pol_vec),
	 M = matconcat([M; Polroots_norm_vec(pol_vec[i])]);
       );
  /* print(matsize(M)); */
  return(M);
};


/* evaluation of precision */
Precision_eval(pol_field, equation, {version=1}) = {
  my(pol_e, ind, base, M, Minv, Binv, pol_mink, Bpol, n_gp, A, prec, field_dim,\
     r, n, m, C, V, n_t);

  field_dim = poldegree(pol_field);
  n = field_dim;
  r = MySignature(pol_field);
  m = (r[1]>0)+1;
  C = max((n-1)*(n+2)\4, 1);
   
  v = Vec(liftall(equation));
  v = vector(length(v), i, vecmax(abs(Vec(v[i]))));
  v = round(log(norml2(vecmax(v))));
    
  my(s = getabstime());   
  if (r[1]==0,
       
      /* M = Complex_mink(pol_field, vecmax(v)+field_dim, r[1], version, field_dim\2-r[1]); */
      M = Complex_mink(pol_field, vecmax(v), r[1], version, field_dim\4);
      ,
      
      M = Complex_mink(pol_field, vecmax(v), r[1], version, (field_dim-r[1])-2*((field_dim-r[1])\4));
      );  
  /* print("Mink computed in ", strtime(getabstime()-s)); */

  default(realprecision, round(vecmax(v)));
  my(s = getabstime());
  pol_mink = Pol_mink_embed(equation, M[2]);
  /* print("pol embed computed in ", strtime(getabstime()-s)); */
  
  V = Polroots_bound_vec(pol_mink);
  /* print("pol roots bound computed in: ", strtime(getabstime()-s)); */
  if (version==1,
    
      my(s = getabstime());
      Minv = M[2]^(-1);
      /* Minv = conj(transpose(M[2])); */
      /* default(realprecision, 10); */
      /* print (Minv * M); */
      Binv = norml2(Minv);
      Bpol = norml2(V);
      n_gp = Binv*Bpol;
      prec = n * ((2*(n-1)+(m-1))*log(2) + log(n_gp*(2+n + 4*C^2) + 1)) / (log(2)*2) - (n*log(n/(2 * Pi* exp(1))) + 2*(n-1)*log(C))/(2*log(2)); ,

      W = vector(length(V), i, V[i]^2/normlp(M[2][i,], 2)^2);
      /* W = vector(length(V), i, V[i]^2/(n * vecmax(abs(M[2][i,]))^2)); */
      n_gp = vecmax(W, &ind);
      
      prec = field_dim*round((log(n_gp))/(2*log(2)));
      );
  return([round(prec), poldegree(equation), round(v)]);
};


/* evaluation of precision : version with the order basis given */
/* Precision_eval_order(pol_field, equation, order_basis, {version=1}) = { */
/*   my(pol_e, ind, base, Morder, M, Minv, Binv, pol_mink, Bpol, n_gp, A, prec,\ */
/*      field_dim, r, n, m, C, V, n_t); */

/*   field_dim = poldegree(pol_field); */
/*   n = field_dim; */
/*   r = MySignature(pol_field); */
/*   m = (r[1]>0)+1; */
/*   C = max((n-1)*(n+2)\4, 1); */
   
/*   v = Vec(liftall(equation)); */
/*   v = vector(length(v), i, vecmax(abs(Vec(v[i])))); */
/*   v = round(log(norml2(vecmax(v)))); */
  
/*   my(s = getabstime());    */
/*   if (r[1]==0, */
       
/*       /\* M = Complex_mink(pol_field, vecmax(v)+field_dim, r[1], version, field_dim\2-r[1]); *\/ */
/*       M = Complex_mink(pol_field, vecmax(v), r[1], version, field_dim\4); */
/*       , */
       
/*       M = Complex_mink(pol_field, vecmax(v), r[1], version, (field_dim-r[1])-2*((field_dim-r[1])\4)); */
/*       /\* M = Complex_mink(pol_field, vecmax(v), r[1], version,  0); *\/ */
       
/*       ); */
/*   Morder = Order_matrix(pol_field, order_basis); */
/*   M = M * Morder; */
/*   /\* print("Mink computed in ", strtime(getabstime()-s)); *\/ */

/*   default(realprecision, round(vecmax(v))); */
/*   my(s = getabstime()); */
/*   pol_mink = Pol_mink_embed(equation, M[2]); */
/*   /\* print("pol embed computed in ", strtime(getabstime()-s)); *\/ */
    
/*   V = Polroots_bound_vec(pol_mink); */
/*   /\* print("pol roots bound computed in: ", strtime(getabstime()-s)); *\/ */
    
/*   if (version==1, */
/*       my(s = getabstime()); */
/*       Minv = M[2]^(-1); */
/*       Binv = norml2(Minv); */
/*       Bpol = norml2(V); */
/*       n_gp = Binv*Bpol; */
/*       prec = n * ((2*(n-1)+(m-1))*log(2) + log(n_gp*(2+n + 4*C^2) + 1)) / (log(2)*2) - (n*log(n/(2 * Pi* exp(1))) + 2*(n-1)*log(C))/(2*log(2)); , */

/*       W = vector(length(V), i, V[i]^2/normlp(M[2][i,], 2)^2); */
/*       n_gp = vecmax(W, &ind); */
/*       n_t = V[ind]^2; */
/*       prec = field_dim*round((log(n_gp))/(2*log(2))); */
/*       ); */
/*   return([round(prec), poldegree(equation), round(v)]); */
/* }; */




/* evaluation of norm */
Norm_eval(pol_field, equation, {version=1})= {
  my(pol_e, ind, base, M, Minv, Binv, pol_mink, Bpol, n_gp, A, prec, field_dim, r, n, m, C, V, n_t);

  field_dim = poldegree(pol_field);
  n = field_dim;
  r = MySignature(pol_field);
  m = (r[1]>0)+1;
  C = max((n-1)*(n+2)\4, 1);
  v = Vec(liftall(equation));
  v = vector(length(v), i, vecmax(abs(Vec(v[i]))));
  v = round(log(norml2(vecmax(v))));
  if (r[1]==0,
      M = Complex_mink(pol_field, vecmax(v)+field_dim\2, r[1], version, field_dim-r[1]);,
      M = Complex_mink(pol_field, vecmax(v)+field_dim\2, r[1], version,  (field_dim-r[1]));
      );
  my(s = getabstime());
  pol_mink = Pol_mink_embed(equation, M[2]);
  V = Polroots_bound_vec(pol_mink);
  if (version==1,
      my(s = getabstime());
      Minv = M[2]^(-1);
      Binv = norml2(Minv);
      Bpol = norml2(V);
      n_gp = Binv*Bpol;  ,
       
      W = vector(length(V), i, V[i]^2/normlp(M[2][i,], 2)^2);
      n_gp = vecmax(W, &ind);
      );
  return(n_gp);
};


/* FUNCTIONS FOR RELATIVE EXTENSIONS */
Rel_pol_mink_bounds(pol, pol_vec, Bg, Me) =  {
  local(m, R);
  m = Rel_polynomial_mink(liftall(pol), pol_vec, Bg, Me);
  R = matrix(poldegree(pol), length(m));
  for(i=1, length(m),
	R[,i] = polrootsbound(m[i]);
      );
  return(R);
};


/* evaluation of precision */
Rel_precision_eval(pol_vec, equation, {version=1}) = {
  my(dg, da, de, dp, pol_e, ind, base, M, Minv, Binv, pol_mink, Bpol, n_gp, A,\
     prec, field_dim, r1, n, m, C, B_vec);
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
  da = dg*de;
  dp = poldegree(equation);
  n = dg;
  C = max((n-1)*(n+2)\4, 1);

  v = Multipol_to_vec(liftall(equation), 3, [dp+1, dg, de]);
  v = vector(length(v), i, norml2(v[i]));
  v = round(log(norml2(vecmax(v)))/2);
  
  my(s = getabstime());

  [M, B_vec, r1] = Abs_complex_minkowski(pol_vec, (v + dg) , version);
  m = (r1>0)+1;
  /* print("mink computed in: ", strtime(getabstime()-s)); */
  
  default(realprecision, v);
  pol_mink = Abs_polynomial_mink(equation, pol_vec, B_vec, version);
  /* print("embeddings and mink computed in: ", strtime(getabstime()-s)); */

  /* print("bounds computed in: ", strtime(getabstime()-s)); */
  V = Polroots_bound_vec(pol_mink);
  default(realprecision, 20);

  if (version==1,
      /* default(realprecision, 50); */
      default(realprecision, v\poldegree(equation));
      M = 1.*M;
      Minv = M^(-1);
      /* print("inv mink computed in: ", strtime(getabstime()-s)); */
      Binv = norml2(Minv);
      Bpol = norml2(V);
      n_gp = Binv*Bpol;
      prec = round(n * ((2*(n-1)+(m-1))*log(2) + log(n_gp*(2+n + 4*C^2) + 1)) / (log(2)*2) - (n*log(n/(2 * Pi* exp(1))) + 2*(n-1)*log(C))/(2*log(2))); ,
      
      W = vector(length(V), i, V[i]/normlp(M[i,], 2));
      n_gp = vecmax(W, &ind);
      
      /* print(V[ind]^2, "\t", normlp(M[ind,])^2); */
      prec = n*round((log(n_gp))/(log(2)));
      );
  return([prec, round(v\3)]);
};


/* /\* /!\ TODO  : MODIFY THIS FUNCTION  *\/ */
/* /\* evaluation of precision *\/ */
/* Rel_precision_eval_order(pol_vec, equation, order_basis, {version=1}) = { */
/*   my(dg, da, de, dp, pol_e, ind, base, M, Minv, Binv, pol_mink, Bpol, n_gp, A,\ */
/*      prec, field_dim, r1, n, m, C, B_vec, Morder); */
/*   dg = poldegree(pol_vec[1]); */
/*   de = poldegree(pol_vec[2]); */
/*   da = dg*de; */
/*   dp = poldegree(equation); */
/*   n = dg; */
/*   C = max((n-1)*(n+2)\4, 1); */

/*   v = Multipol_to_vec(liftall(equation), 3, [dp+1, dg, de]); */
/*   v = vector(length(v), i, norml2(v[i])); */
/*   v = round(log(norml2(vecmax(v)))); */

/*   my(s = getabstime()); */

/*   [M, B_vec, r1] = Abs_complex_minkowski(pol_vec, v\3 + dg, version); */
/*   m = (r1>0)+1; */
/*   /\* print("mink computed in: ", strtime(getabstime()-s)); *\/ */
/*   Morder = Order_matrix(polvec[1], order_basis); */
  
  
/*   default(realprecision, v\2); */
/*   pol_mink = Abs_polynomial_mink(equation, pol_vec, B_vec, version); */
/*   /\* print("embeddings and mink computed in: ", strtime(getabstime()-s)); *\/ */

/*   V = Polroots_bound_vec(pol_mink); */
/*   /\* print("bounds computed in: ", strtime(getabstime()-s)); *\/ */
    
/*   if (version==1, */
/*       Minv = M^(-1); */
/*       /\* print("inv mink computed in: ", strtime(getabstime()-s)); *\/ */
/*       Binv = norml2(Minv); */
/*       Bpol = norml2(V); */
/*       n_gp = Binv*Bpol; */
/*       prec = round(n * ((2*(n-1)+(m-1))*log(2) + log(n_gp*(2+n + 4*C^2) + 1)) / (log(2)*2) - (n*log(n/(2 * Pi* exp(1))) + 2*(n-1)*log(C))/(2*log(2))); , */
     
/*       W = vector(length(V), i, V[i]^2/normlp(M[i,], 2)^2); */
/*       n_gp = vecmax(W, &ind); */
/*       prec = n*round((log(n_gp))/(2*log(2))); */
/*       ); */
/*   return([prec, round(v\3)]); */
/* }; */
