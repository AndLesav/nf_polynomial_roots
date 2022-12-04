MyGramMatrix(M)=
  {
    return(mattranspose(M)*M);
  };

/* naive babai */
MyGSO_naive(M, r, prec)={
  my(G = MyGramMatrix(M), B, Q, R, mu, MU, N, dim, nrows, ncols);
  [nrows, ncols] = matsize(M);
  B = matrix(nrows, ncols);
  N = matrix(1, ncols);
  B[,1] = M[,1];
  default(realprecision, prec);
  for (i = 2, ncols,
	 B[,i] = M[,i];
	 for (j = 1, i-1,
		mu = (B[,j]~)*B[,i];
	      B[,i] -= (1.0*mu/norml2(B[,j]))*B[,j];
	      );
       );
  return(B);
};


/* using recursive formulas and not computing everything */
MyGSO(M, r, prec)={
  my(B, Q, R, mu, MU, N, dim, nrows, ncols);
  [nrows, ncols] = matsize(M);
  MU = matrix(ncols, ncols);
  N = vector(ncols);
  B = matrix(r, ncols);
  /* initialising with first column */
  B[,1] = M[1..r, 1];
  N[1] = norml2(M[,1]);
  default(realprecision, prec);
  for (i = 2, ncols,
	 B[,i] = M[1..r,i];
       N[i] = norml2(M[,i]);
       for (j = 1, i-1,
	      /* determining mu_{i,j} with rec. formula  */
	      mu = (M[,i]~)*M[,j];
	    mu -= sum(k=1, j-1, MU[k,j]*MU[k,i]*N[k]);
	    N[i] -= ((1.*mu^2)/N[j]);
	    mu /= (1.)*N[j];
	    MU[j, i] = mu;
	    B[,i] -= mu*B[,j];
	    /* N[i] -= mu^2*N[j]; */
	    );
       );
  return([B, matdiagonal(N)*MU, N]);
};

/* using recursive formulas and not computing everything following
   L2 paper by Nguyen and Stehlé */
MyGSO_L2(M, r, prec)={
  my(B, Q, R, mu, MU, N, dim, nrows, ncols);
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
       /* N[i] = norml2(M[,i]); */
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
  return([B, R, N]);
};


/* babai using QR rep - assume matrix has been completed  */
MyBabai(v, M, QR)=
  {
   my(w, r);
   r = matsize(M)[2];
   w = v;
   for (i=0, r-1,
	  /* print(i+1); */
	  w = w - round(((w~)*QR[1][,r-i])/QR[2][r-i,r-i])*M[,r-i];
	/* print(round(log(norml2(w[2..length(w)])+1))); */
	);
   return(v-w);
  };


/* babai using QR rep - assume matrix has been completed  */
MyDecodeBabai(embeddings, M, QR, prec)=
  {
    my(w, r, field_dim, S, C, Y, c, s, e, q, E, Q, b);

    S = matsize(M);
    
    field_dim = S[2];
   
    default(realprecision, prec+field_dim^2);
   
    C = ((field_dim-1)*(field_dim+2))\4;
    /* C = 1; */
    r = S[1]-S[2];
    /* print(r); */
    Y = vector(r, i, round(embeddings[i]<<prec));
    /* default(realprecision, 10); */
    /* print(norml2(Y)); */
    /* default(realprecision, 2*prec + 10*field_dim^2); */
    v = vector(field_dim);
    v = concat(Y, v);
    w = v~;
    /* print(w); */
    
    /* default(realprecision, prec + field_dim^2); */
    s = vector(field_dim);
    /* E = vector(field_dim); */
    /* Q = vector(field_dim); */
    
    
    /* first inner product : multiplication by one or two scalars */
    /* d = 0; */
    /* for (i =1, r, */
    /* 	   d += Y[i]*QR[1][i,field_dim]; */
    /* 	 ); */
    /* d=d/QR[2][field_dim,field_dim]; */
    /* q = round(d); */
    /* c = d; */

    c = (Y*QR[1][1..r,field_dim])/QR[2][field_dim,field_dim];
    q = round(c);

    e = abs(c-q);
    /* print(e); */
    b = (e < field_dim\10);
    /* b=1; */
    
    my(i = 0);
    while (b && (i < field_dim-1),
	   /* print(i); */
	   default(realprecision, prec+field_dim^2);
	   w = w - q*M[,field_dim-i];
	   
	   i++;

	   c = ((w~)*QR[1][,field_dim-i])/QR[2][field_dim-i,field_dim-i];
	   
	   q = round(c);
	   e += abs(c-q);

	   default(realprecision, 10);
  	   
	   b = (e < field_dim\10);
	   
	   /* print(e); */
	   /* i++; */
	   
	   );
	  
    if (b, w = w - q*M[,1];);
    /* if (e>100, print(e)); */
    return([b, w[r+1..field_dim+r]/C]);

   /*  /\* s = 0; *\/ */
   /*  my(i = 1); */
   /*  for (i=1, field_dim-1, */
   /* 	   default(realprecision, prec + field_dim^2); */
   /* 	 /\* print(i+1); *\/ */
   /* 	 c = ((w~)*QR[1][,field_dim-i])/QR[2][field_dim-i,field_dim-i]; */
   /* 	 q = round(c); */
	
   /* 	/\* print(c); *\/ */
   /* 	/\* print(q); *\/ */
   /* 	e = c-q; */
   /* 	w = w - q*M[,field_dim-i]; */
   /* 	s[field_dim-i] = e*QR[2][field_dim-i,field_dim-i]; */
   /* 	E[field_dim-i] = e; */
   /* 	Q[field_dim-i] = q; */
   /* 	/\* default(realprecision, 10); *\/ */
   /* 	/\* print(e); *\/ */
   /* 	/\* s[field_dim-i] = (w~)*QR[1][,field_dim-i]; *\/ */
   /* 	default(realprecision, 10); */
   /* 	/\* print(e); *\/ */
   /* 	if(vecsum(abs(E)) > 0.25, print("error too big at step: ", i+1); i = field_dim;); */
   /* 	/\* if (i==0, print(w)); *\/ */
   /* 	/\* print(s[field_dim-i]); *\/ */
   /* 	); */
   /* /\* print(w[1]); *\/ */
   /* default(realprecision, 10); */
   /* print("quotient: ", log((norml2(w)/norml2(s)))); */
   /* print("norm in GSO: ", round(log(norml2(s)))); */
   /* /\* print("norm coeff :", norml2(E)); *\/ */
   /* print("norm coeff :", vecsum(abs(E))); */
   /* /\* print(norml2(Q)); *\/ */
   /* /\* print(round(log(norml2(w[r+1..field_dim+r]/C)))); *\/ */
   /* return(w[r+1..field_dim+r]/C); */
  };



/* babai using QR rep - assume matrix has been completed  */
MyDecodeBabai_1(embeddings, M, QR, prec, {heur=1})=
  {
    my(w, r, field_dim, S, C, Y, c, s, e, q, E, Q, b, err, W);

    S = matsize(M);
    
    field_dim = S[2];
   
    default(realprecision, prec+field_dim);
   
    C = ((field_dim-1)*(field_dim+2))\4;
    r = S[1]-S[2];
    Y = vector(r, i, round(embeddings[i]<<prec));
    v = vector(field_dim);
    v = concat(Y, v);
    w = v~;
    
    /* default(realprecision, prec + field_dim^2); */
    s = vector(field_dim);
        
    W = vector(field_dim, i, Y*QR[1][1..r,i]);

    c = W[field_dim]/QR[2][field_dim,field_dim];
    q = round(c, &err);
    /* print(c); */
    b = 1;

    if (heur, 
	e = abs(c-q);
	b = (e < 1/4);
	);
    /* b=1; */
    
    my(i = 0);
    while (b && (i < field_dim-1),

	   default(realprecision, prec+field_dim);
	   W -= Vec(q*QR[2][1..field_dim-i-1,field_dim-i], field_dim);
	   
   	   w -= q*M[,field_dim-i];
	   
	   i++;
	   
	   c = W[field_dim-i]/QR[2][field_dim-i, field_dim-i];

	   q = round(c, &err);
   
	   if (heur, 
	       
	       default(realprecision, 10);
	       b = (abs(c-q) < 1/4);
	       );
	   
	   /* if (!b, print("CONDITIONS IN BABAI NOT VERIFIED");); */
	   
	   );

    if (b, w = w - q*M[,1];);

    print("***********************************");
    
    return([b, w[r+1..field_dim+r]/C]);
    
  };




/* babai using GSO rep   */
MyDecodeBabai_2(embeddings, M, G, prec, {heur=1})=
  {
    my(w, r, field_dim, S, C, Y, c, s, e, q, E, Q, b, err, W);

    S = matsize(M);
    
    field_dim = S[2];
   
    default(realprecision, prec+field_dim);
   
    C = ((field_dim-1)*(field_dim+2))\4;
    r = S[1]-S[2];
    Y = vector(r, i, round(embeddings[i]<<prec));
    v = vector(field_dim);
    v = concat(Y, v);
    w = v~;

    s = vector(field_dim);
        
    W = vector(field_dim, i, Y*G[1][,i]);
    
    c = W[field_dim]/G[3][field_dim];
    q = round(c);
    
    b = 1;

    if (heur, 
	e = abs(c-q);
	b = (e < 1/4);
	);
    
    my(i = 0);
    while (b && (i < field_dim-1),
	   default(realprecision, prec);
	   W -= Vec(vector(field_dim-i-1, j, q*G[2][j,field_dim-i]), field_dim);
   	   w -= q*M[,field_dim-i];
	   
	   i++;
	   
	   c = W[field_dim-i]/G[3][field_dim-i];
	   q = round(c);

	   if (heur, 
	       
	       default(realprecision, 10);
	       b = (abs(c-q) < 1/4);
	       );
	   	   
	   );

    if (b, w = w - q*M[,1];);
    return([b, w[r+1..field_dim+r]/C]);
  };



/* /\* ind = sum e_i*dv[i]^i<--> [e_i]_i  *\/ */
/* IndexToExponent_multi(ind, dv)= */
/*   { */
/*    my(i, j, e, d, n, p); */
/*    n = length(dv); */
   
/*    j = ind-1; */
/*    e = vector(n); */
/*    /\* print(e) *\/ */
/*    for(i=1, n, */
/* 	 d = dv[i+1..n]; */
/*        p = vecprod(d); */
/*        [e[i], j] = divrem(j, p)+1; */
/*        ); */
/*    return(e); */
/*   }; */


/* /\* ind = sum e_i*dv[i]^i<--> [e_i]_i  *\/ */
/* ExponentToIndex_multi(e, dv)= */
/*   { */
/*    my(i, j, ind, d, n, p); */
/*    n = length(dv); */
/*    ind = 0; */
/*    /\* print(e) *\/ */
/*    for(i=1, n, */
/*        	 d = dv[i+1..n]; */
/*        p = vecprod(d); */
/*        ind += (e[i]-1)*p; */
/*        /\* [e[i], j] = divrem(j, p)+1; *\/ */
/*        ); */
/*    return(ind+1); */
/*   }; */


/* FirstSameInd(e, old_e)= */
/*   { */
/*    my(b=1, k=0); */
   
/*    while(b && (k < length(e)), */
/* 	 k += 1; */
/* 	 /\* print(k, " ", length(old_e), " ", length(e)); *\/ */
/* 	 b=prod(j=1, length(old_e[k]), e[k]!=old_e[k][j]);	  */
/* 	 ); */
/*    return([b, k]); */
/*   }; */



/* FirstIndBefore_multi(e, k, dv)= */
/*   { */
/*    my(b=0, l=-1); */
/*    while(!b && (l < (k-1)), */
/* 	 l += 1; */
/* 	 /\* if(p==45,print(l);); *\/ */
/* 	 b = (e[k-l]!=dv[k-l]); */
/* 	 ); */
/*    return([b, k-l]); */
/*   }; */


/* NextInd_multi(e, dv, old_e)= */
/*   { */
/*    my(e_new, k, l, b, b1, k1, n); */
/*    n = length(e); */
/*    /\* q = (length(e)-1)/2; *\/ */
/*    /\* b1 = prod(j=1, q, e[2*j-1+r1]==e[2*j+r1]); *\/ */
/*    e_new = e; */
/*    if (prod(j=1, length(e), e[j]==dv[j]), return(e);); */
   
/*    [b, k] = FirstSameInd(e, old_e); */
/*    /\* print([b, k]); *\/ */
/*    if (b, return(e)); */
   
/*    while (!b, */
/* 	  /\* if(e_new[1]==45,print(e_new " " k);); *\/ */
/* 	  [b1, k1] = FirstIndBefore_multi(e_new, k, dv); */
/* 	  /\* print([b1, k1]); *\/ */
/* 	  if(!b1, */
/* 	     return(e_new); , */
	     
/* 	     e_new[k1]=e_new[k1]+1; */
/* 	     for(l=k1+1, n, e_new[l]=1); */
	     
/* 	     ); */
/* 	  [b, k] = FirstSameInd(e_new, old_e); */
/* 	  /\* print([b, k]); *\/ */
/* 	  ); */
   
/*    return(e_new); */
/*   }; */


/* UpdatePair_conj(e, d) = */
/*   { */
/*    /\* print(e); *\/ */
/*    /\* (e[1]==d && e[2]==d, return([0, [1, 1]]); , *\/ */
/*    if(e[1]==e[2], return([0, [e[1], e[2]]]); , */
/*       e[1]<e[2], return([1, [e[1]+1, e[1]+1]]); , */
/*       e[2]<e[1], return([1, [e[1], e[1]]]);  */
/*       ); */
/*   }; */


/* UpdatePair_conj_cond(e, d) = */
/*   { */
/*    /\* print(e); *\/ */
/*    if (e[1]==d && e[2]==d, return([0, [1, 1]]); , */
/*        e[1]==e[2], return([1, [e[1]+1, e[2]+1]]); , */
/*       e[1]<e[2], return([1, [e[1]+1, e[1]+1]]); , */
/*       e[2]<e[1], return([1, [e[1], e[1]]]);  */
/*       ); */
/*   }; */



/* NextConj_multi_simple(e, dv, r1) = */
/*   { */
/*    my(n, r2, b, a, e_new, e_pair); */
/*    n = length(e); */
/*    r2 = (n-r1)/2; */
   
/*    if (prod(j=1, n, e[j]==dv[j]), return(e);); */
 
/*    e_new = e; */
   
/*    for (i=1, r2, */
	  
/* 	  [a, e_pair] = UpdatePair_conj(e_new[r1+2*i-1..r1+2*i], dv[r1+2*i-1]); */

/* 	/\* print(a); *\/ */
/* 	if (a, */

/* 	    e_new[r1+2*i-1] = e_pair[1]; */
/* 	    e_new[r1+2*i] = e_pair[2]; */

/* 	    for (j = i+1, r2, */
/* 		   e_new[r1+2*j-1] = 1; */
/* 		 e_new[r1+2*j] = 1; */
/* 		 ); */
	    
/* 	    return(e_new); */
/* 	    ); */

/* 	); */

   
/*    /\* for (i=0, r2-1, *\/ */
	  
/*    /\* 	  [a, e_pair] = UpdatePair_conj_cond(e_new[n-2*i-1..n-2*i], dv[n-2*i-1]); *\/ */
/*    /\* 	/\\* print(e_pair); *\\/ *\/ */
/*    /\* 	e_new[n-2*i-1] = e_pair[1]; *\/ */
/*    /\* 	e_new[n-2*i] = e_pair[2]; *\/ */
	
/*    /\* 	if (a, *\/ */
	    
/*    /\* 	    return(e_new)); *\/ */
/*    /\* 	); *\/ */
/*    /\* print(e_new); *\/ */
/*    /\* print(dv); *\/ */

/*    /\* for (i=0, r1-1, *\/ */

/*    /\* 	  if (e_new[r1-i]!=dv[r1-i], *\/ */
/*    /\* 	      return(e_new);, *\/ */

/*    /\* 	      e_new[r1-i]=1; *\/ */
	      
/*    /\* 	      ); *\/ */
	  
/*    /\* 	); *\/ */

/*    return(e); */
   
/*   }; */



Real_basis_multi(pol_field, prec, prec_add, ne) =
  {
   my(field_dim, a, basis, M, R, Rabs);

   field_dim = poldegree(pol_field);
   
   default(realprecision, prec + field_dim^2 + prec_add);
   
   /* my(t = getabstime()); */
   R = polrootsreal(pol_field);
   /* print(R); */
   Rabs = vecsort(R, abs, 4);
   
   M = matrix(ne, field_dim);
     
   for (i = 1, ne,
   	  a = Rabs[i];
   	basis = vector(field_dim, i, a^(i-1));
   	M[i,] = basis;
   	);
   /* a = Rabs[length(Rabs)]; */
   /* basis = vector(field_dim, i, a^(i-1)); */
   /* M[ne,] = basis; */
   return(M);
  };


/* Create base lattice [B|C * Id]~ embedded in RR */
/* use several embeddings ; assume ne < r_1 */
Real_basis_lattice_multi(pol_field, prec, Me, uni, {spec_lll=1})=
  {
   my(field_dim, a, basis, B, M, uni_new, C);
   field_dim = poldegree(pol_field);
   /* C = 2^(prec \ (2) ); */
   C = max((field_dim-1)*(field_dim+2)\4, 1);
   /* C = 1; */
   default(realprecision, prec\3 + field_dim);
   basis = Me[1,];
   /* print(basis); */

   basis = Me[1,];
   /* for (i = 2, matsize(Me)[1], */
   /* 	  for (j = 1, matsize(Me)[2],  */
   /* 		 basis[j] *= Me[i,j]; */
   /* 	       ); */
   /* 	); */

   B = vector(field_dim, i, -round(basis[i]<<prec));

   for (i = 2, matsize(Me)[1],
	  /* print(B); */
	  basis = Me[i,];
	B = matconcat([B ; vector(field_dim, j, -round(basis[j]<<prec))]);
	);

   /* print("Basis creation took: " strtime(getabstime()-t)); */
   /* default(realprecision, 20); */
   
   M = B*uni;
   M = matconcat([M; C*uni]);
 
   my(s = getabstime());
   if (spec_lll==1,
       print("computing with spec lll");
       M = MySpecLLL(M); ,
              
       M = qflll(M, 3);
       );
   print("LLL computed in: ", strtime(getabstime()-s));
   
   uni_new = M[matsize(Me)[1]+1..field_dim+matsize(Me)[1],]/C;
   /* uni_new = M[2..field_dim+1,]/C; */
   
   return([M, uni_new]);
};


Test_solve_equation_real_multi(pol_field, equation, prec, Me, lattice, prec_add, {heur=1, ns=poldegree(equation)})=
  {
   local(q, r, g1, gemb, pol, field_dim, R, S, M, P, pol1, ne, dv);
   field_dim = poldegree(pol_field);
   ne = matsize(Me)[1];
   old_e = vector(ne, j, []);

   default(realprecision, round(prec\5)+2*field_dim+prec_add\10);
   my(s = getabstime());
   /* QR = matqr(matsupplement(lattice)); */
   QR = MyGSO_L2(lattice, ne, prec+2*field_dim);
   print("QR computation: ",  strtime(getabstime()-s));

   default(realprecision, round(prec\2+prec_add+field_dim));
   my(s = getabstime());
   pol = Pol_mink_embed_nf(equation, Me); /* put equation in RR */
   print("embedding of equation: ",  strtime(getabstime()-s));

   [M, P] = FF_basis(pol_field, 32, field_dim\4);
   pol1 = FF_pol_embedding_fam_nf(equation, M, P); /* put eq in FF */
   
   default(realprecision, round(prec\3)+field_dim+prec_add\5);
 
   R = vector(length(pol));
   dv = vector(length(pol));
   
   for(i=1, length(pol),
	 R[i] = polrootsreal(pol[i]);
       dv[i] = length(R[i]);
       /* print(R[,i]); */
       );
   print(dv);
   print("roots of equations: ",  strtime(getabstime()-s));
   S = [];
   
   for (i = 1, vecprod(dv),
	  /* if (i%10==0, print(i)); */
	  /* print(length(S)); */
	e = IndexToExponent_multi(i, dv);
	e = NextInd_multi(e, dv, old_e);
	i = ExponentToIndex_multi(e, dv);
	r = vector(ne, j, R[j][e[j]]); 
	my(s = getabstime());
	[b, g1] = MyDecodeBabai_2(r, lattice, QR, prec, heur);
	g1 = g1~;

	/* print("Decoding done in: ", strtime(getabstime()-s)); */
	if (b, 
	    my(s = getabstime());
	    
	    gemb = FF_eval_pol(g1, pol1, M, P);
                    
	    g1 = Mod(Polrev(g1, y), pol_field);
       
	    my(s = getabstime());
	    
	    if(FF_test_vec(gemb, P),
	       /* print("TRUE"); */
	       /* if(r==0,  */
	       /* print("verification done in: ", strtime(getabstime()-s)); */
	       S = (concat(S, [g1]));
	       if(length(S)==ns && ns!=0, return(S));
	       e_new = vector(length(e), j, 1);
	       e_new[1] = e[1]+1;
	       i = ExponentToIndex_multi(e_new, dv)-1;
	       for(j=1, ne,
		     old_e[j] = matconcat([old_e[j], e[j]])[1,];
		   );
	       
	       i = ExponentToIndex_multi(e, dv); ,
	       
	       /* print("FALSE"); */

	       );
	    );
	if (number_solutions != 0,
	    if (length(S)==number_solutions, return(S););
	    );
	/* print(length(S)); */
	);
   /* print(S, g "\n"); */
   return(S);
  };
