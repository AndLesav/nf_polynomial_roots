/* ######################### FUNCTIONS FOR LATTICES ######################### */

/*########################################################################### */
/*                         GSO + REDUCTION + DECODING                         */
/*########################################################################### */

/* specific lll for power-basis */
MySpecLLL(B, {index=2}) =  {
  my(r, C, M, field_dim = matsize(B)[2]);
   
  C = max((field_dim-1)*(field_dim+2)\4, 1);
  r = matsize(B)[1]-matsize(B)[2];
  M = Mat(B[1..(r+index), 1..index]);

  if (index > 1,
      M[index+r, index] = C*round(norml2(M[,index-1])^2);
      );

  M = qflll(M, 3);

  for (i=index+1, field_dim,
	 M[i+r-1,i-1] = C;
       M = matconcat([M, B[1..i+r,i]]);
       for (j=1, i-2, M[,i] += (B[1..i+r,j+1]/C)*M[j+r,i-1];);
       if (i < field_dim, M[i+r, i] = C*round(norml2(M[,i-1])^2););
       my(t = getabstime());
       M = qflll(M, 3);
       );
   
  return(M);
};


/* returns potential element given an approx. embedding in RR or CC 
   uses Kannan embedding
 */
Test_decode(embeddings, prec, lattice) = {
    local(S, x, y, h, i, v, field_dim, M, C, r);
    S = matsize(lattice);
    field_dim = S[2];
    C = ((field_dim-1)*(field_dim+2))\4;
    r = S[1]-S[2];
    Y = vector(r, i, round(embeddings[i]<<prec));
    v = vector(field_dim+1);
    v = concat(Y, v);
    v[field_dim+r+1] = 2^(prec\2);
    v = v~;
    M = matconcat([lattice, v]);

    my(s=getabstime());
    M = qflll(M, 3);
    /* print("Time for decoding: ", strtime(getabstime()-s)); */

    M = M[,field_dim+1]~;
    return(M[r+1..field_dim+r]/C);
  };


/* using recursive formulas and not computing everything following
   L2 paper by Nguyen and Stehlé --
   more efficient than matqr from gp (?)*/
MyGSO_L2(M, r, prec) = {
  my(B, Q, Q1, R, mu, MU, N, dim, nrows, ncols);
  [nrows, ncols] = matsize(M);
  MU = matrix(ncols, ncols);
  R = matrix(ncols, ncols);
  N = vector(ncols);
  B = matrix(r, ncols);

  /* initialising with first column */
  B[,1] = M[1..r, 1];
  N[1] = norml2(M[,1]);
  R[1,1] = N[1];
  MU[1,1] = 1;
  
  /* starting the other columns */
  default(realprecision, prec);
  for (i = 2, ncols,
	 B[,i] = M[1..r,i];
       for (j = 1, i-1,
	      /* determining mu_{i,j} with rec. formula  */
	      mu = (M[,i]~)*M[,j];
	    mu -= sum(k=1, j-1, MU[k,j]*R[k,i]);
	    R[j,i] = mu;
	    
	    mu /= (1.0*R[j,j]);
	    MU[j, i] = mu;
	    B[,i] -= mu*B[,j];
	    );
       mu = (M[,i]~)*M[,i];
       mu -= sum(k=1, i-1, MU[k,i]*R[k,i]);
       R[i,i] = mu;
       N[i] = R[i,i];
       );
  Q1 = vector(length(B), i, B[,i]/N[i]);
  return([B, R, N, Q1]);
  /* return([B, R, N]); */
};



/* babai using GSO rep 
   embeddings : approx embedding of element 
   M : basis lattice 
   G : vector describing GSO as returned by MyGSO_L2
   prec : precisionat which computations are done
   heur : 0 for certified method and 1 for heuristic method
*/
MyDecodeBabai_2(embeddings, M, G, prec, {heur=1}) = {
  my(w, r, field_dim, S, C, Y, c, s, e, q, E, Q, b, err, W);
  S = matsize(M);
  field_dim = S[2];

  default(realprecision, prec\3+field_dim\2);
  C = ((field_dim-1)*(field_dim+2))\4;
  r = S[1]-S[2];
  Y = vector(r, i, round(embeddings[i]<<prec));
  v = vector(field_dim);
  v = concat(Y, v);
  w = v~;
  Q = vector(field_dim);
  
  /* W = vector(field_dim, i, Y*G[1][,i]); */
  /* c = W[field_dim]/G[3][field_dim]; */
   
  /* c = Y*G[4][field_dim]; */
  /* q = round(c); */

  my(s = getabstime());
  c = Y*G[1][,field_dim]/G[3][field_dim];
  q = round(c);


  Q[field_dim]=q;
  
  b = 1;

  if (heur,    /* test if error is not too large */
      e = abs(c-q);
      b = (e < 1/16);
      );
  
  my(i = 0);
  while (b && (i < field_dim-1),
	 default(realprecision, prec\2);

	 w -= q*M[,field_dim-i];

	 /* W -= Vec(vector(field_dim-i-1, j, q*G[2][j,field_dim-i]), field_dim); */

	 i++;

	 /* c = W[field_dim-i]/G[3][field_dim-i]; */
	 
	 
	 c = Y*G[1][,field_dim-i];
	 for (j=0, i-1,
		c -= Q[field_dim-j]*G[2][field_dim-i, field_dim-j];
	      );
	 c /= G[3][field_dim-i];
	 
	   
	 q = round(c);
	 Q[field_dim-i]=q;

	 if (heur, /* test if error is not too large */
	     default(realprecision, 10);
	     b = (abs(c-q) < 1/16);
	     );
	 );

  if (b, w = w - q*M[,1];);
  /* print("stop after ", i+1, " steps"); */
  return([b, w[r+1..field_dim+r]/C]);
};



MyDecodeBabai_nbsteps(embeddings, M, G, prec, {heur=1}) = {
  my(w, r, field_dim, S, C, Y, c, s, e, q, E, Q, b, err, W);
  S = matsize(M);
  field_dim = S[2];

  default(realprecision, prec\2+field_dim\2);
  C = ((field_dim-1)*(field_dim+2))\4;
  r = S[1]-S[2];
  Y = vector(r, i, round(embeddings[i]<<prec));
  v = vector(field_dim);
  v = concat(Y, v);
  w = v~;
  Q = vector(field_dim);
  
  /* W = vector(field_dim, i, Y*G[1][,i]); */
  /* c = W[field_dim]/G[3][field_dim]; */
   
  /* c = Y*G[4][field_dim]; */
  /* q = round(c); */

  my(s = getabstime());
  c = Y*G[1][,field_dim]/G[3][field_dim];
  q = round(c);


  Q[field_dim]=q;
  
  b = 1;

  if (heur,    /* test if error is not too large */
      e = abs(c-q);
      b = (e < 1/8);
      );
  
  my(i = 0);
  while (b && (i < field_dim-1),
	 default(realprecision, prec\2);

	 w -= q*M[,field_dim-i];

	 /* W -= Vec(vector(field_dim-i-1, j, q*G[2][j,field_dim-i]), field_dim); */

	 i++;

	 /* c = W[field_dim-i]/G[3][field_dim-i]; */
	 
	 
	 c = Y*G[1][,field_dim-i];
	 for (j=0, i-1,
		c -= Q[field_dim-j]*G[2][field_dim-i, field_dim-j];
	      );
	 c /= G[3][field_dim-i];
	 
	   
	 q = round(c);
	 Q[field_dim-i]=q;

	 if (heur, /* test if error is not too large */
	     default(realprecision, 10);
	     b = (abs(c-q) < 1/8);
	     );
	 );

  if (b, w = w - q*M[,1];);
  /* print("stop after ", i+1, " steps"); */
  return([b, w[r+1..field_dim+r]/C, i+1]);
};




/* ################################################################################ */
/*                  FUNCTIONS FOR THE ABSOLUTE CASE                                 */
/* ################################################################################ */


/* Create base lattice [B|C * Id]~ embedded in RR */
Create_basis_lattice(pol_field, prec) = {
  my(field_dim, a, basis, B, M, uni_new, C, R, Rabs, ind);
  field_dim = poldegree(pol_field);
  C = max((field_dim-1)*(field_dim+2)\4, 1);
   
  default(realprecision, prec\2 + field_dim^2 + prec);
   
  my(t = getabstime());
  R = polrootsreal(pol_field);
  Rabs = abs(R);
   
  a = vecmax(Rabs, &ind);
  a = R[ind];
  basis = vector(field_dim, i, a^(i-1));
  kill(a);
    
  B = vector(field_dim, i, -round(basis[i]<<prec));
  default(realprecision, prec \ 2.5 +1);
   
  my(t = getabstime());
  M = matconcat([B; C*matid(field_dim)]);

  return(M);
};



/* Create base lattice [B|C * Id]~ embedded in CC */
Create_basis_lattice_complex(pol_field, prec) = {
  my(field_dim, a, basis, B, B_real, B_im, M, uni_new, C, R, Rabs, ind);
  field_dim = poldegree(pol_field);
  C = max((field_dim-1)*(field_dim+2)\4, 1);
   
  default(realprecision, prec\2 + field_dim^2 + prec);
  my(t = getabstime());
  R = polroots(pol_field);
  Rabs = abs(R);
  a = vecmax(Rabs, &ind);
  a = R[ind];
  basis = vector(field_dim, i, a^(i-1));
  kill(a);

  B = vector(field_dim, i, -basis[i]<<prec);
  B_real = vector(field_dim, i, round(real(B[i])));
  B_im = vector(field_dim, i, round(imag(B[i])));
  M = matconcat([B_real; B_im]);
  M = matconcat([M; C*matid(field_dim)]);
  return(M);
};


/* Creat basis lattice [B|C * Id]~ embedded in RR and REDUCE IT  */
Real_basis_lattice(pol_field, prec, uni, prec_add, {spec_lll=1, order_basis=[]}) = {
  my(field_dim, a, basis, B, M, uni_new, C, R, Rabs, ind, basis_eq);
  field_dim = poldegree(pol_field);
  C = max((field_dim-1)*(field_dim+2)\4, 1);

  default(realprecision, prec + field_dim^2 + prec_add\2);
  my(t = getabstime());
  R = polrootsreal(pol_field);
  Rabs = abs(R);
  /* print("roots computation took: " strtime(getabstime()-t)); */
    
  a = vecmax(Rabs, &ind);	/* max seems better experimentally */
  a = R[ind];
  basis_eq = vector(field_dim, i, a^(i-1)); 

  if (length(order_basis) != 0,
      basis = vector(length(order_basis), i,
		     Element_embedding(Vecrev(order_basis[i]), basis_eq)); ,
      basis = basis_eq
      );
  
  kill(a);
  kill(basis_eq);
  
  B = vector(field_dim, i, -round(basis[i]<<prec));
  default(realprecision, prec \ 2.5 +1);
  
  my(t = getabstime());
  M = B*uni;			/* pre-reduce */
  M = matconcat([M; C*uni]);
  kill(B);
  my(t = getabstime());
  if (spec_lll==1,
      M = MySpecLLL(M); ,
      M = qflll(M, 3);
      );
  uni_new = M[2..field_dim+1,]/C;
  /* print("LLL took: " strtime(getabstime()-t)); */
  return([M, uni_new, basis]);
};



/* Create base lattice [B | C*Id]~ embedded in CC */
Complex_basis_lattice(pol, prec, uni, prec_add, {spec_lll=1, order_basis=[]}) = {
  my(field_dim, a, basis, B, M, uni_new,  ind, b, basis1);
  field_dim = poldegree(pol);
  C = (field_dim-1)*(field_dim+2)\4;

    default(realprecision, round(prec\2)+field_dim+prec_add);
  my(t = getabstime());
  R = polroots(pol);
  /* print("Time taken for roots computation: ", strtime(getabstime()-t)); */
  Rabs = vector(length(R), i, abs(R[i]));
  a = vecmax(Rabs, &ind);
  a = R[ind];
  basis_eq = vector(field_dim, i, a^(i-1)); 

  if (length(order_basis) != 0,
      basis = vector(length(order_basis), i,
		     Element_embedding(Vecrev(order_basis[i]), basis_eq)); ,
      basis = basis_eq
      );

  kill(a);
  kill(basis_eq);
  
  default(realprecision, round(prec/3)+field_dim*field_dim+prec_add);
  B = vector(field_dim, i, -basis[i]<<prec);
  B_real = vector(field_dim, i, round(real(B[i])));
  B_im = vector(field_dim, i, round(imag(B[i])));
  /* print("Basis creation took: " strtime(getabstime()-t)); */
   
  kill(a);

  M = matconcat([B_real; B_im]);
  M = M*uni;			/* pre-reduce with previous unitary transform */
  M = matconcat([M; C*uni]);

  my(t = getabstime());
  if (spec_lll,
      M = MySpecLLL(M); ,
      M = qflll(M, 3);
      );
  uni_new = M[3..field_dim+2,]/C;
  return([M, uni_new, basis]);
};


/* initialisation */
Real_basis_lattice_init(pol, prec_add, {order_basis=[]})=
  {
    local(uni);
    uni = matid(poldegree(pol));
    return(Real_basis_lattice(pol, 1, uni, prec_add, 0, order_basis));
  };


/* initialisation */
Complex_basis_lattice_init(pol, prec_add, {order_basis=[]})=
  {
    local(uni);
    uni = matid(poldegree(pol));
    return(Complex_basis_lattice(pol, 1, uni, prec_add, 0, order_basis));
  };

  
/* goes from precision "precision" to precision "precision*q^k"  with step q */
Real_basis_lattice_update_geom(pol, prec, prec_add, q, k, uni) = {
  my(lattice, basis);
  for(i=1, k,
	print (i, "over", k);
	prec = ceil(prec*q);
      [lattice, uni, basis] = Real_basis_lattice(pol, prec, uni, prec_add, 0);
      );
  return([lattice, uni, basis, prec]);
};


/* goes from precision "precision" to precision "precision*q^k"  with step q */
Complex_basis_lattice_update_geom(pol, prec, prec_add, q, k, uni) = {
  my(lattice, basis);
  for(i=1, k,
	prec = ceil(prec*q);
      [lattice, uni, basis] = Complex_basis_lattice(pol, prec, uni, prec_add, 0);
      );
  return([lattice, uni, basis, prec]);
};

/* compute new reduced lattice [B | C*Id] in RR from scratch */
New_basis(pol, target_prec, prec_add, {geom=0, spec_lll=1, order_basis=[]}) = {
  my(prec, B, uni, prec_s);
  [B, uni, prec] = Real_basis_lattice_init(pol, prec_add, order_basis);
  s = getabstime();
  if (geom,
      prec_s=1;
      expo = floor(log(target_prec/prec_s)/log(2));
      [B, uni, basis, prec] = Real_basis_lattice_update_geom(pol, prec_s, prec_add, 2, expo, uni);
      [B, uni, basis] = Real_basis_lattice(pol, target_prec, uni, prec_add, 0);
      ,
      
      [B, uni, basis] = Real_basis_lattice(pol, target_prec, matid(poldegree(pol)), prec_add, spec_lll, order_basis);
      );   
  return([B, uni, basis, target_prec]);
};

/* compute new reduced lattice [B | C*Id] in CC from scratch */
New_basis_complex(pol, target_prec, prec_add, {geom=0, order_basis=[]}) = {
  my(prec, B, uni, basis);
  [B, uni, basis] = Complex_basis_lattice_init(pol, prec_add, order_basis);

  prec = 1;

  if (geom,
      expo = floor(log(target_prec)/log(2));
      [B, uni, basis, prec] = Complex_basis_lattice_update_geom(pol, prec, prec_add,\
								2, expo, uni);
      [B, uni, basis] = Complex_basis_lattice(pol, target_prec, uni, prec_add); ,
       
      [B, uni, basis] = Complex_basis_lattice(pol, target_prec, uni, prec_add, 1, order_basis);
      );   
  return([B, uni, basis, target_prec]);
};




/* ################################################################################ */
/*                  FUNCTIONS FOR THE RELATIVE CASE                                 */
/* ################################################################################ */

/* Create base lattice [Bg | C*Id]~ embedded in CC */
Rel_complex_basis_lattice(pol_vec, Bg, precision, uni, {spec_lll=0}) = {
  my(field_dim, a, basis, B, M, uni_new, ind, C);
  field_dim = poldegree(pol_vec[1]);
  C = max((field_dim-1)*(field_dim+2)\4, 1);
  default(realprecision, round(precision/3)+field_dim*field_dim);
   
  my(t = getabstime());
  basis = Bg;
  B = vector(field_dim, i, -basis[i]<<precision);
  B_real = vector(field_dim, i, round(real(B[i])));
  B_im = vector(field_dim, i, round(imag(B[i])));
  /* print("Basis creation took: " strtime(getabstime()-t)); */

  M = matconcat([B_real; B_im]);
  M = M*uni;
  M = matconcat([M; C*uni]);
   
  my(t = getabstime());
  default(realprecision, round(precision));

  if (spec_lll,
      M = MySpecLLL(M); ,
      M = qflll(M, 3);
      );
   
  uni_new = M[3..field_dim+2,]/C;
  return([M, uni_new, basis]);
};


Rel_partial_decode(mink, Me_inv, lattice, precision, {ind=0, gind=[]})={
  my(emb, g, n, test);
  emb = Me_inv*mink~;
  g = vector(length(emb));
  n = vector(length(emb));
  for(i=1, length(emb),
	if(i!=ind, 
	   g[i] = Test_decode([real(emb[i]), imag(emb[i])], precision, lattice);,
	   g[i] = gind;
	   );
      default(realprecision, 5);
      );
  return(g)
};


Rel_partial_decode_qr(mink, Me_inv, lattice, QR, precision, {heur=1, ind=0, gind=[]})={
  my(emb, g, n, test, h, i);

  g = vector(length(mink));

  i = 0;

  bt = 1;
  
  while (bt && i < length(g),
	 i += 1;
	 if(i!=ind,
	    /* emb = Me_inv[i,]*mink~; */
	    emb = sum(j=1, length(mink), mink[j][i]);
	    
	    [bt, h] = MyDecodeBabai_2([real(emb), imag(emb)], lattice, QR, \
				      precision, heur);

	    g[i] = h~; ,
	    
	    g[i] = gind~;
	    );
	 );
  return([bt,g])
};


/* test one coord. with QR  */
Rel_test_one_coord_qr(mink, Me_inv, lattice, QR, precision,  {heur = 1, ind = 1}) = {
  my(emb, gt, n, b, quo, a);
  default(realprecision, precision\3);
  emb = sum(j=1, length(mink), mink[j][ind]);
  /* emb = Me_inv[ind,]*mink~; */
  
  [b, gt] = MyDecodeBabai_2([real(emb), imag(emb)], lattice, QR, precision, heur);
  return([b, b, gt]);
};
