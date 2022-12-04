/* FUNCTIONS DEALING WITH EMBEDDINGS IN FINITE FIELDS */


/* FIRST : GENERAL FUNCTIONS / ABSOLUTE CASE  */

/* Several basis of K with finite fields */
FF_basis(pol_field, s, N) = {
  local(B, M, a, p, j, R, n, P);
  n = poldegree(pol_field);
  M = matrix(N, n);
  P = [];
  for (i = 1, N,
	 R = [];
       while (length(R)==0, 
	      p = randomprime(2^s);
	      R = polrootsmod(pol_field, p);
	      );
       P = concat(P, p);
       j = 1+random(length(R));
       a = R[j];
       B = vector(n, i, a^(i-1));
       M[i,] = B;
       );
  return([M, P]);
};


/* element embedding in finite fields */
FF_element_embedding(g, B, p) = {
  return(B*(Mod(g, p)~));
};


/* element embedding - several embeddings */
FF_element_embedding_fam(g, M, P) = {
  local(h, G);
  G = [];
  for (i = 1, length(P),
	 h = Mod(g, P[i]);
       G = concat(G, M[i,]*(h~));
       );
  return(G);
};


/* pol is pol in x with coeff. pol in y */
FF_pol_embedding(pol, B, p)= {
  local(v, w, k, n, m, V);
  w = Vec(pol);
  m = length(B);
  v = vector(poldegree(pol)+1, k, FF_element_embedding(Vecrev(w[k], m), B, p));
  return(Pol(v));
};


/* pol is in x with coeff in nf as pol in y mod p */
FF_pol_embedding_nf(pol, B, p) = {
  local(q);
  q = lift(pol);		/* lift coeff as pol in y */
  return(FF_pol_embedding(q, B, p));
};



/* pol is pol in x with coeff. pol in y */
FF_pol_embedding_fam(pol, M, P) = {
  local(V);
  V = [];
  for (i = 1, length(P),
	 V = concat(V, FF_pol_embedding(pol, M[i,], P[i]));
       );
  return(V);
};


/* pol is in x with coeff in nf as pol in y mod p */
FF_pol_embedding_fam_nf(pol, M, P) = {
  local(q);
  q = lift(pol);		/* lift coeff as pol in y */
  return(FF_pol_embedding_fam(q, M, P));
};


FF_eval_pol(g, pol_fam, M, P) = {
  my(h, R);
  R = [];
  for (i=1, length(P),
	 h = FF_element_embedding(g, M[i,], P[i]);
       R = concat(R, subst(pol_fam[i], x, h));
       );
  return(R);
};


FF_test_vec(v, P) = {
  local(b, i);
  b = 1;
  i = 1;
  while (b &&  i < length(P)+1,
	 b *= (v[i]==Mod(0, P[i]));
	 i += 1;
	 );
  return(b);
};
  

/* ****************************************************************************** */
/* ********************  functions for relative ext. L/K  *********************** */
/* ****************************************************************************** */

/* several bases of L/K in finite fields */
FF_rel_basis(pol_vec, s, N) = {
  my(a, b, dg, de, da, Bg, Be, pe, B_vec, R, S, P, p);
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
  da = dg*de;

  B_vec = vector(N);

  P = [];
   
  for (i = 1, N,
	 b = 0;
       while (!b,
	       
	      R = [];

	      while (length(R)==0,

		     p = randomprime(2^s);

		     R = polrootsmod(pol_vec[1], p);
		     );

	      j = 1+random(length(R));
	      a = R[j];

	      Bg = vector(dg, i, a^(i-1));

	      pe = FF_pol_embedding_nf(pol_vec[2], Bg, p);
	      S = polrootsmod(pe, p);
	      b = (length(S)!=0);

	      );
	
       P = concat(P, p);
       j = 1+random(length(S));
       a = S[j];
       Be = powers(a, de-1);
       B_vec[i] = [Bg, Be];
       );
      
  return([B_vec, P]);
};

/* embedding of x in several FF when seen as vector in Kg^de */
FF_rel_element_embedding(g, pol_vec, B, p) = {
  local(k, v, w);
	
  w = vector(length(g), j, FF_element_embedding(g[j], B[1], p));
	   
  return(B[2]*w~);
   
};


/* embedding of x in several FF when seen as vector in Kg^de */
FF_rel_element_embedding_fam(g, pol_vec, B_vec, P) = {
  local(k, v, w, p);

  v = vector(length(P));

  for (i=1, length(P),

	 p = P[i];
	
       w = vector(length(g), j, FF_element_embedding(g[j], B_vec[i][1], p));
	
       v[i] = B_vec[i][2]*w~;
	
       );
   
  return(v);
   
};


/* pol is pol. in x with coeff. pol in (QQ[y])[z] ?? */
FF_rel_pol_embedding(pol, pol_vec, B, p) = {
  local(w, de, dg, dp);
   
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
  dp = poldegree(pol);
   
  w = Multipol_to_vec(pol, 3, [dp+1, de, dg]);
  w = vector(dp+1, i, FF_rel_element_embedding(w[i], pol_vec, B, p));
   
  return(Polrev(w, variable(pol)));
};

/* pol is pol. in x with coeff. pol in (QQ[y])[z] ?? */
FF_rel_pol_embedding_nf(pol, pol_vec, B, p) = {
  return(FF_rel_pol_embedding(liftall(pol), pol_vec, B, p));
};


/* pol is pol in x with coeff. pol in y */
FF_rel_pol_embedding_fam(pol, pol_vec, B_vec, P) = {
  local(V);

  V = [];

  for (i = 1, length(P),
	 V = concat(V, FF_rel_pol_embedding(pol, pol_vec, B_vec[i], P[i]));
       );

  return(V);
};

/* pol is in x with coeff in nf as pol in y mod p */
FF_rel_pol_embedding_fam_nf(pol, pol_vec, B_vec, P) = {
  local(q);
   
  q = liftall(pol);		/* lift coeff as pol in y */

  return(FF_rel_pol_embedding_fam(q, pol_vec, B_vec, P));
};



FF_rel_eval_pol_fam(g, pol_vec, pol_fam, B_vec, P) =
  {
    my(h, R);

    R = [];
   
    for (i=1, length(P),
	  
	   h = FF_rel_element_embedding(g, pol_vec, B_vec[i], P[i]);

	 R = concat(R, subst(pol_fam[i], x, h));
	  
	 );
   
    return(R);
  };

