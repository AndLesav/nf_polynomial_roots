/*  Generate random HNF of the form [C1 | ... | Cr | C * I_n]~
    Ci is a column vector with random entries in [0,Det]
*/
HNF_gen(n, Det, C,{r=1}) = {
   my(M = C*matid(n));
   for (i = 1, r, 
	  v = vector(n, i, random(Det));
	M = matconcat([v ;M]);
	);
   return (M);
};


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

shuffle_indexes(n) = {
  my (I, L);
  I = List(vector(#l, i, i));
  L = List(vector(#l));
  for (i = 1, #l,
	 a = random(#I) + 1;
       L[i] = I[a];
       listpop(~I, a);
       );
  return(L);
};

list_shuffle(l) = {
  my (a, I, L);
  I = List(vector(#l, i, i));
  L = List(vector(#l));
  for (i = 1, #l,
	 a = random(#I) + 1;
       L[i] = l[I[a]];
       listpop(~I, a);
       );
  return(L);
};

list_shuffle_pair(l, g) = {
  my (a, I, L, G);
  I = List(vector(#l, i, i));
  L = List(vector(#l));
  G = List(vector(#l));
  for (i = 1, #l,
	 a = random(#I) + 1;
       L[i] = l[I[a]];
       G[i] = g[I[a]];
       listpop(~I, a);
       );
  return([L, G]);
};


Vec_shuffle(V, r1_g, r1) = {
  my(W = V, L, G);
  for (i = 1, r1,
	 W[i] = list_shuffle(V[i]);
       );
  for (i = 1, (#V - r1)\2,
	 [L, G] = list_shuffle_pair(V[r1 + 2*i-1], V[r1 + 2*i]);
       W[r1+2*i-1] = L;
       W[r1+2*i] = G;
       );
  return(W);
}


/* ************************************************************************** */
/*                       FUNCTIONS CHANGE REPRESENTATION                      */
/* ************************************************************************** */

/* goes from basis of order given as rational polynomials to matrix */
Order_matrix(p, order_basis) = {
  my (M, dim);
  dim = poldegree(p);
  M = [];
  for (i = 1, dim,
	 M = matconcat([M, Vecrev(order_basis[i], dim)~]);
       );
  return (M);
}


/* ************************************************************************** */
/*                FUNCTIONS CHANGE REPRESENTATION RELATIVE EXT                */
/* ************************************************************************** */

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

myProdvec(v) = {
  my (n = length(v), w = vector(length(v)+1));
  w[1] = 1;
  for (i=1, length(v),
	 w[i+1] = w[i]*v[n-i+1];
       );
  return(w);
}
  

/* ind = sum e_i*dv[i]^i<--> [e_i]_i  */
IndexToExponent_multi(ind, dv, {pv =0}) =  {
  my(i, j, e, d, n, p);
  n = length(dv);
  if (pv ==0, p = myProdvec(dv);, p = pv);
  j = ind-1;
  e = vector(n);
  /* print(e) */
  for(i=1, n,
      [e[i], j] = divrem(j, p[n-i+1])+1;
      );
  return(e);
};


/* ind = sum e_i*dv[i]^i<--> [e_i]_i  */
ExponentToIndex_multi(e, dv, {pv= 0}) = {
  my(i, j, ind, d, n, p);
  n = length(dv);
  if (pv ==0, p = myProdvec(dv);, p = pv);
  ind = 0;
  return(vecsum(vector(n, i, (e[i]-1)*p[n-i+1]))+1)
  /* for(i=1, n, */
  /*     /\* 	d = dv[i+1..n]; *\/ */
  /*     /\* p=vecprod(d); *\/ */
  /*     ind += (e[i]-1)*p[n-i+1]; */
  /*     ); */
  /* return(ind+1); */
};

NextExp_multi(~e, dv) = {
  /* my (e_new = e); */
  my (i = length(e));
  while (i > 0 && e[i] == dv[i],
	 e[i] = 1;
	 i -= 1;
	 );
  if (i > 0, e[i] += 1 );
  return (e);
};

FirstSameInd(e, old_e) = {
  my(b=1, k=0);
   
  while(b && (k < length(e)),
	k += 1;
	j = 0;
	b=prod(j=1, length(old_e[k]), e[k]!=old_e[k][j]);
	);
  return([b, k]);
};


FirstIndBefore_multi(e, k, dv) = {
  my(b=1, l=-1);
  /* b = prod(l=1, k-1, e[k-l]!=dv[k-l]); */
  while(b && (l < (k-1)),
	l += 1;
	b = (e[k-l]==dv[k-l]);
	);
  return([!b, k-l]);
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


TestExpoConj(e, p, n, r1) = {
  my(q, b);
  q = (length(e)-r1)/2;
  b = prod(j=1, q, e[2*j-1+r1]==e[2*j+r1]);
  return(b);
};



BatchMatMult(M, R) = {
  my(W, V = vector(length(R), i, []));
  for (i=1, length(R),
       W = vector(length(R[i]), j, M[,i]*R[i][j];);
       V[i] = List(W);
       );
  return (V);
};






/* ******************************************************************************** */
/*               functions for good / bad field in Pari/Gp                          */
/* ******************************************************************************** */


get_maxf(dim) = {
  my(maxf = 1);
  if (dim >= 45, maxf = 32; ,
      dim >= 30, maxf = 16; ,
      dim >= 15, maxf = 8;
      );
  return (maxf);
};

get_nbprimes(pol_field, {deg_eq = 1}) = {
  my(B, nb);
  if (deg_eq != 1, B = poldegree(pol_field) * deg_eq;,
      B = poldegree(pol_field) * 20;
      );
  if (B <= 128, nb = 5; ,
      B <= 1024, nb = 20; ,
      B <= 2048, nb = 65;,
      nb = 100;
      );
  return (max(nb, 2*poldegree(pol_field)));
};


max_inertial_degree(pol_field, p) = {
  my(F);
  F = factormod(pol_field, p, 1);
  return (vecmax(F));
};

is_good_prime(pol_field, p, maxf) = {
  my(S, id = max_inertial_degree(pol_field, p));
  if (id==poldegree(pol_field),
      S = [1, id];, 
      S = [id <= maxf, id];
      );
  return(S);
};
  
find_best_inertial_degree(pol_field, {deg_eq=1}) = {
  my(disc, dim, maxf, best_f, init, p, S);
  disc = abs(poldisc(pol_field));
  dim = poldegree(pol_field);
  maxf = get_maxf(dim);
  best_f = 1;
  init = 0;
  p = 1;
  while (!init,
	 p = nextprime(p+1);
	 if ((disc % p) != 0,
	     S = is_good_prime(pol_field, p, maxf);
	     init = S[1];
	     best_f = S[2];
	     );
	 );
  my(count = get_nbprimes(pol_field, deg_eq)-1);
  while ((count > 0) && (best_f != dim),
	 p = nextprime(p+1);
	 if ((disc % p) != 0,
	     count--;
	     S = is_good_prime(pol_field, p, maxf);
	     /* print(S); */
	     if (S[1],
		 best_f = max(S[2], best_f);
		 );
	     );
	 );
  return (best_f);
};

is_good_field(pol_field, {deg_eq = 1}) = {
  my (dim = poldegree(pol_field), best_f, bool);
  best_f = find_best_inertial_degree(pol_field, deg_eq);
  return ([best_f == dim, best_f]);
};
 
