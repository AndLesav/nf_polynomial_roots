/* FUNCTIONS RELATED TO EMBEDDINGS IN RR or CC */


/* FUNCTIONS IN THE ABSOLUTE SITUATION */

/* Create base embedded in RR */
Real_basis(pol_field, prec, index) = {
  my(field_dim, a, basis, B, M, uni_new, C);
  field_dim = poldegree(pol_field);
  default(realprecision, round(prec/3)+field_dim);
  my(t = getabstime());
  a = polrootsreal(pol_field)[index];
  basis = vector(field_dim, i, a^(i-1));
  kill(a);
  return(basis);
};


/* real roots + "real mink" in matrix form */
Real_mink(pol_field, r1, prec) = {
  my(field_dim, R, M);
  field_dim = poldegree(pol_field);
  default(realprecision, round(prec/3)+field_dim);
  R = polrootsreal(pol_field);
  M = matrix(r1, field_dim);
  for(i=1, length(R),
	M[i,] = powers(R[i], field_dim-1);
      );
  return([R, M]);
};

 
/* complex roots + complex mink in matrix form */
Complex_mink(pol_field, prec, r1, {version, ne=poldegree(pol_field)-r1}) = {
  my(field_dim, R, M);
  field_dim = poldegree(pol_field);
   
  default(realprecision, round(prec/3)+field_dim);
  my(s = getabstime());
  if (r1!=0 && ne==0,
      R = polrootsreal(pol_field); ,
      R = polroots(pol_field);
      );
  /* print("roots computed in ", strtime(getabstime()-s)); */
   
  if (version,
      M = matrix(field_dim);
      for(i=1, field_dim,
	    M[i,] = powers(R[i], field_dim-1);
	  ); ,

      M = matrix(r1+ne\2, field_dim);
      for(i=1, r1,
	    M[i,] = powers(R[i], field_dim-1);
	  );
      for(i=1, ne\2,
	    M[r1+i,] = powers(R[r1+2*i], field_dim-1);
	  );
      );
  return([R, M]);
};


/* inverse mink by interpolation */
Complex_mink_inv(pol_field, R) = {
  my(field_dim, M, p);
  field_dim = poldegree(pol_field);
  M = matid(field_dim);
  for(i=1, field_dim,
	v = vector(field_dim)~; 
      v[i] = 1;
      p = Vecrev(polinterpolate(R, v))~;
      M[,i] = p;
      );
  return(M);
};



/* embedding into RR / CC of element is in vector rep  */
Element_embedding(g, basis) = {
  local(k);
  return(sum(k=1, length(g), basis[k]*g[k]));
};

  
/* embedding into RR / CC pol is pol in x with coeff. pol in y */
Pol_embedding(pol, basis) = {
  local(v, w, k);
   
  w = Vec(pol);
  v = vector(poldegree(pol)+1, k, 
	     Element_embedding(Vecrev(w[k], length(basis)), basis));
   
  return(Pol(v));
};

/*  image of pol. under mink of pol is pol in x with coeff. pol in y */
Pol_mink_embed(pol, M) = {
  local(v, w, k, pol_vec);
   
  pol_vec = vector(matsize(M)[1], i, Pol_embedding(lift(pol), M[i,]));
  
  return(pol_vec);
};


/* pol is in x with coeff in nf as pol in y mod p */
Pol_embedding_nf(pol, basis) = {
  local(q);
    
  q = lift(pol);		/* lift coeff as pol in y */
    
  return(Pol_embedding(q, basis));
};


/*  pol is pol in x with coeff. pol in y */
Pol_mink_embed_nf(pol, M) = { 
  local(v, w, k, pol_vec);

  pol_vec = vector(matsize(M)[1]);
   
  for(i = 1, length(pol_vec),
	pol_vec[i] = Pol_embedding_nf(pol, M[i,]);
      );
    
  return(pol_vec);
};



/* FUNCTIONS FOR THE RELATIVE SITUATION */

/* create basis of [ground, ext] in CC */
Complex_basis_creation_rel(pol_vec, precision, index_vec) = {
    my(a, b, dg, de, Bg, Be, pe);
    dg = poldegree(pol_vec[1]);
    de = poldegree(pol_vec[2]);
   
    a = polroots(pol_vec[1])[index_vec[1]];
    Bg = powers(a, dg-1);

    pe = Pol_embedding_nf(pol_vec[2], Bg);
    a = polroots(pe)[index_vec[2]];
    Be = powers(a, de-1);

    return([Bg, Be]);
};


/* relative mink Ke/Kg by basis + matrix representation */
Rel_complex_minkowski(pol_vec, precision) = {
  my(a, b, dg, de, Bg, Be, pe, r1, M, r1_g);
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
   
  default(realprecision, round(precision));
   
  a = polrootsreal(pol_vec[1]);

  if (length(a)!=0, 
      b = vector(length(a), i, abs(a[i]));
      b = vecmax(b, &ind);
      a = a[ind];
      r1_g = 1;,

      a = polroots(pol_vec[1]);
      b = vector(length(a), i, abs(a[i]));
      b = vecmax(b, &ind);
      a = a[ind];
      r1_g=0;
      );
   
  Bg = powers(a, dg-1);
   
  pe = Pol_embedding_nf(pol_vec[2], Bg);
  R = polroots(pe);
  M = matid(de);
  for(i=1, length(R),
	M[i,] = powers(R[i], de-1);
      );

  if (imag(Bg[2])!=0,
      r1 = 0; ,
       
      r1 = sum(i=1,length(R),imag(R[i])==0);
      );
       
  return([Bg, M, r1, r1_g]);
};


/* relative mink Ke/Kg by basis + matrix representation */
Abs_complex_minkowski(pol_vec, precision, {version=1}) = {
  my(a, b, dg, de, da, Bg, Be, pe, r1, M, N, v, w, B_vec, c, Me);
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
  da = dg*de;
  r1 = 0;
   
  default(realprecision, round(precision));
   
  B_vec = vector(dg);
  c=1;
  a = polroots(pol_vec[1]);

  Bg = powers(a[1], dg-1);
  pe = Pol_embedding_nf(pol_vec[2], Bg);

  R = polroots(pe);

  if (imag(Bg[2])!=0,
      r1 += 0; ,
       
      r1 += sum(i=1,length(R),imag(R[i])==0);
      );
   
   
  M = matrix(de, da);
  Me = matid(de);
  for(i=1, length(R),
	v = powers(R[i], de-1);
      Me[i,] = v;
      w = [];
      
      for (j=1, length(v),
	     w = concat(w, Bg*v[j]);
	   );
      M[i,] = w;
      );
  B_vec[c] = [Bg, Me];
  c += 1;
    
  N = M;
    
  if (version==1,
      le=length(a); ,
      le=length(a)\2+1;
      );
  for (k=2, le,
	 Me = matid(de);

       Bg = powers(a[k], dg-1);
      
       pe = Pol_embedding_nf(pol_vec[2], Bg);
       R = polroots(pe);

       if (imag(Bg[2])!=0,
	   r1 += 0; ,
	    
	   r1 += sum(i=1,length(R),imag(R[i])==0);
	   );

	
       M = matrix(de, da);
   
       for(i=1, length(R),
	     v = powers(R[i], de-1);
	   Me[i,] = v;
	   w = [];
	   for (j=1, length(v),
		  w = concat(w, Bg*v[j]);
		);
	   M[i,] = w;
	   );
       B_vec[c] =  [Bg, Me];
       /* print(B_vec); */
       c +=1;
       N = matconcat([N;M]);
	
       );	  
  return([N, B_vec, r1]);
};



/* embedding of g viewed as an elt of Ke/Kg i.e. vector in Kg^de */
Rel_element_embedding(g, pol_vec, basis_vec) = {
  local(k, v);
  v = vector(length(g), i, Element_embedding(g[i], basis_vec[1]));
   
  return(sum(k=1, length(g), basis_vec[2][k]*v[k]));
};


/* rel. mink. of x viewed as an elt of Ke/Kg i.e. vector in Kg^de */
Rel_element_mink(g, pol_vec, Bg, Me) = {
  local(k, v);
  v = vector(length(g), i, Element_embedding(g[i], Bg));
   
  return(Me*(v~));
};


/* pol is pol. in x with coeff. pol in (QQ[y])[z] ?? */
Rel_polynomial_embedding(pol, pol_vec, basis_vec) = {
  local(w, de, dg, dp);
   
  dg = poldegree(pol_vec[1]);
  de = poldegree(pol_vec[2]);
  dp = poldegree(pol);
   
  w = Multipol_to_vec(pol, 3, [dp+1, de, dg]);
  w = vector(dp+1, i, Rel_element_embedding(w[i], pol_vec, basis_vec));
   
  return(Polrev(w, variable(pol)));
};


/* rel. mink. of pol in Ke/Kg[x]  */
Rel_polynomial_mink(pol, pol_vec, Bg, Me) = {
    local(m, w, de, dg, dp);
   
    dg = poldegree(pol_vec[1]);
    de = poldegree(pol_vec[2]);
    dp = poldegree(pol);
   
    w = Multipol_to_vec(pol, 3, [dp+1, de, dg]);
    m = Rel_element_mink(w[1], pol_vec, Bg, Me);
    for(i=2, length(w),
	  m = matconcat([m, Rel_element_mink(w[i], pol_vec, Bg, Me)]);
	);
    w = vector(matsize(m)[1], i, Polrev(m[i,], variable(pol)));
   
    return(w);
};


Abs_polynomial_mink(pol, pol_vec, B_mink, {version=1}) = {
  local(v, l);
  v = [];
  if (version,
      l = length(B_mink);,
      l = length(B_mink)\2+1;,
      );
  for (i=1, l,
	 /* print(Rel_polynomial_mink(liftall(pol), pol_vec, B_mink[i][1], B_mink[i][2])); */
	 my(s = getabstime());
       v = concat(v, Rel_polynomial_mink(liftall(pol), pol_vec, B_mink[i][1], B_mink[i][2]));
       /* print("one embedding computed in: ", strtime(getabstime()-s)); */
       );
  return(v);
};


/* roots of mink. of pol. */
Rel_pol_mink_roots(pol, pol_vec, Bg, Me, r1, prec) =  {
  local(m, R, dg);

  dg = poldegree(pol_vec[1]);
  my(s = getabstime());
  default(realprecision, prec + dg);
  m = Rel_polynomial_mink(liftall(pol), pol_vec, Bg, Me);
  /* printf("time taken for mink: %s\n", strtime(getabstime()-s)); */
  R = matrix(poldegree(pol), length(m));
  R = vector(length(m), i, []);
    
  default(realprecision, prec);
   
  for(i=1, r1,
	R[i] = List(polrootsreal(real(m[i])));
      );
   
  for(i=r1+1, length(m),
	s = getabstime();
      R[i] = List(polroots(m[i]));
      /* printf("time taken for roots comp is: %s\n", strtime(getabstime()-s)); */
      
      );

  return(R);
};
