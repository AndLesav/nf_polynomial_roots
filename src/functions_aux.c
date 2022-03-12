/* signature */
MySignature(pol)= {
  my(r, r1, r2);
  R = polrootsreal(pol);
  r1 = length(R);
  r2 = (poldegree(pol)-r1)\2;
  return([r1, r2]);
};

/* update equation */
Update_equation(pol_field, equation, S)={
  my(eq=equation);
  for (i=1, length(S),
	 eq /= (x-S[i]);
       );
  return(eq);
};

/* median value */
MyMedian(v)={
  my(n = length(v));
  if (n % 2 == 0, 
      return((v[n\2] + v[n\2+1])/2); ,
      return(v[n\2+1]);
      );
};


/* Square sum */
MySquareSum(v, c)={
  my(n = length(v));
  return(vecsum(vector(n, i, (v[i]-c)^2)));
};


/* least square to sort in two groups */
MyVecSort_two_groups(v)={
  my(n = length(v));
  my(V = vector(n-1));
  my(w = vecsort(v));
  my(c, m1, m2, w1, w2, ind, M1, M2, m);

  M1 = vector(n-1);
  M2 = vector(n-1);
   
  for (i=1, n-1,
	 c = (w[i]+w[i+1])/2;

       w1 = w[1..i];
       w2 = w[i+1..n];
       m1 = MyMedian(w1);
       M1[i] = m1;
       m2 = MyMedian(w2);
       M2[i] = m2;
       V[i] = MySquareSum(w1, m1) + MySquareSum(w2, m2);
	   
       );

  m = vecmin(V, &ind);
  return([M1[ind], M2[ind], w[1..ind], w[ind+1..n] ]);
};

/* ******************************************************************************** */
/* FUNCTIONS CHANGE REPRESENTATION RELATIVE EXT */
/* ******************************************************************************** */

/* goes from pol. rep to vec. rep. */
Multipol_to_vec(pol, nv, deg_vec) = {
  my(v);
  v = Vecrev(pol, deg_vec[1]);
  if(nv==1,
     return(v),
     v = vector(length(v), i, Multipol_to_vec(v[i], nv-1, deg_vec[2..nv]));
     return(v);
     );
};

/* goes from vec. to pol. rep. */
Vec_to_multipol(v, nv, {varg=z, vare=y}) = {
  my(pol);
  pol = vector(length(v), i, Polrev(v[i], varg));
  pol = Polrev(pol, vare);
  return(pol);
};



/* ******************************************************************************** */
/*               functions for search in relative methods                          */
/* ******************************************************************************** */

/* ind = sum e_i*dv[i]^i<--> [e_i]_i  */
IndexToExponent_multi(ind, dv) =  {
  my(i, j, e, d, n, p);
  n = length(dv);
   
  j = ind-1;
  e = vector(n);
  /* print(e) */
  for(i=1, n,
	d = dv[i+1..n];
      p = vecprod(d);
      [e[i], j] = divrem(j, p)+1;
      );
  return(e);
};


/* ind = sum e_i*dv[i]^i<--> [e_i]_i  */
ExponentToIndex_multi(e, dv) = {
  my(i, j, ind, d, n, p);
  n = length(dv);
  ind = 0;
  /* print(e) */
  for(i=1, n,
	d = dv[i+1..n];
      p = vecprod(d);
      ind += (e[i]-1)*p;
      /* [e[i], j] = divrem(j, p)+1; */
      );
  return(ind+1);
};


FirstSameInd(e, old_e) = {
  my(b=1, k=0);
   
  while(b && (k < length(e)),
	k += 1;
	b=prod(j=1, length(old_e[k]), e[k]!=old_e[k][j]);
	);
  return([b, k]);
};


FirstIndBefore_multi(e, k, dv) = {
  my(b=0, l=-1);
  while(!b && (l < (k-1)),
	l += 1;
	b = (e[k-l]!=dv[k-l]);
	);
  return([b, k-l]);
};


NextInd_multi(e, dv, old_e) = {
  my(e_new, k, l, b, b1, k1, n);
  n = length(e);
  e_new = e;
  if (prod(j=1, length(e), e[j]==dv[j]), return(e););
   
  [b, k] = FirstSameInd(e, old_e);
  if (b, return(e));
   
  while (!b,
	 [b1, k1] = FirstIndBefore_multi(e_new, k, dv);
	 if(!b1,
	    return(e_new); ,
	    e_new[k1]=e_new[k1]+1;
	    for(l=k1+1, n, e_new[l]=1);
	    );
	 [b, k] = FirstSameInd(e_new, old_e);
	 );
  return(e_new);
};


UpdatePair_conj(e, d) = {
  if(e[1]==e[2], return([0, [e[1], e[2]]]); ,
     e[1]<e[2], return([1, [e[1]+1, e[1]+1]]); ,
     e[2]<e[1], return([1, [e[1], e[1]]]); 
     );
};


UpdatePair_conj_cond(e, d) = {
  if (e[1]==d && e[2]==d, return([0, [1, 1]]); ,
      e[1]==e[2], return([1, [e[1]+1, e[2]+1]]); ,
      e[1]<e[2], return([1, [e[1]+1, e[1]+1]]); ,
      e[2]<e[1], return([1, [e[1], e[1]]]); 
      );
};




NextConj_multi_simple(e, dv, r1) = {
  my(n, r2, b, a, e_new, e_pair);
  n = length(e);
  r2 = (n-r1)/2;
   
  if (prod(j=1, n, e[j]==dv[j]), return(e););

  e_new = e;
  for (i=1, r2,
	 [a, e_pair] = UpdatePair_conj(e_new[r1+2*i-1..r1+2*i], dv[r1+2*i-1]);
       if (a,
	   e_new[r1+2*i-1] = e_pair[1];
	   e_new[r1+2*i] = e_pair[2];
	   for (j = i+1, r2,
		  e_new[r1+2*j-1] = 1;
		e_new[r1+2*j] = 1;
		);
	   return(e_new);
	   );
       );
  return(e);
};


TestExpoConj(e, p, n, r1, {old_e=vector(n, j, [0])}) = {
  my(q, b);
  q = (length(e)-r1)/2;
  b = prod(j=1, q, e[2*j-1+r1]==e[2*j+r1]);
  return(b);
};
