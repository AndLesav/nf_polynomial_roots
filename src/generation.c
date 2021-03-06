/* CREATION OF POLYNOMIALS */

/* POL. DEFINING FIELDS */

/* generates irreducible pol. with at least one real root */
Pol_field_creation_real(field_dim, size_coeff, {var=y})= {
  my(p, r1, R);
  r1 = 0;
  while(r1==0,
	p = vector(field_dim, i, random(2*2^size_coeff)-2^size_coeff);
	p = Pol(p, var);
	p += var^field_dim;
	while(!polisirreducible(p),
	      p = vector(field_dim, i, random(2*2^size_coeff)-2^size_coeff);
	      p = Pol(p, var);
	      p += var^field_dim;
	      );
	R = polrootsreal(p);
	r1 = length(R);
	);
  return(p);
};


/* generates irreducible pol. with at least one complex root */
Pol_field_creation_CC(field_dim, size_coeff, {var=y})=
  {
    my(p, r1, R);
    r1 = field_dim;
    while(r1==field_dim,
	  p = vector(field_dim, i, random(2*2^size_coeff)-2^size_coeff);
	  p = Pol(p, var);
	  p += var^field_dim;
	  while(!polisirreducible(p),
		p = vector(field_dim, i, random(2*2^size_coeff)-2^size_coeff);
		p = Pol(p, var);
		p += var^field_dim;
		);
	  R = polrootsreal(p);
	  r1 = length(R);
	  );
    return(p);
  };


/* generates irreducible pol. with no real root */
Pol_field_creation_complex(field_dim, size_coeff, {var=y})=
  {
    my(p, r1, R);
    if (Mod(field_dim,2)!=0,
	print("ERROR WRONG DIM");
	return 0; ,
       
	r1 = 1;
	while(r1!=0,
	      p = vector(field_dim, i, random(2*2^size_coeff)-2^size_coeff);
	      p = Pol(p, var);
	      p += var^field_dim;
	      while(!polisirreducible(p),
		    p = vector(field_dim, i, random(2*2^size_coeff)-2^size_coeff);
		    p = Pol(p, var);
		    p += var^field_dim;
		    );
	      R = polrootsreal(p);
	      r1 = length(R);
	      );
	);
    return(p);
  };


/* generates irreducible polynomial */
Pol_field_creation(field_dim, size_coeff, {var=y}) =
  {
    my(p);
    p = sum(k=0, field_dim-1, (random(2*2^size_coeff)-2^size_coeff)*var^k);
    p += var^field_dim;
    while(!polisirreducible(p),
	  p = sum(k=0, field_dim-1, (random(2*2^size_coeff)-2^size_coeff)*var^k);
	  p += var^field_dim;
	  );
    return(p);
  };


/* POL. DEFINING EQUATIONS */

/* P defined by integral coefficients */
Pol_creation_integral(degree, size_coeff, {var=x})={
  my(p);
  p = random((2^size_coeff+1)*var^degree);
  return(p);
};


/* P with coefficient in number field */
Pol_eq_creation_field(degree, size_coeff, pol_field, {var_p=x, var_f=y}) = {
  my(q, coeff);
  q = 0;
  coeff = (2^size_coeff+1)*var_f^poldegree(pol_field);
  for(l=0, degree,
	q += random(coeff)*var_p^l;
      );
  return(q);
};


/* PROBLEM CREATION */
/* create equation of the form h = q(g) with q in ZZ[X] or nf[X] */
Equation_creation(pol_field, size_root, degree_pol, size_pol, type) = {
  my(g, h, q, field_dim);

  field_dim = poldegree(pol_field);
   
  g = vector(field_dim, k, random(2*2^size_root)-2^size_root);
  g = Mod(Polrev(g, y), pol_field);
   
  if (type==0,
      q = Pol_creation_integral(degree_pol, size_pol);,
      q = Pol_eq_creation_field(degree_pol, size_pol, pol_field);
      );
   
  h = subst(q, x, g);
   
  h = Vecrev(lift(h));
   
  return([g, h, q])
};


/* create element in nf[X]  */
Random_el_nf(pol_field, size)={
  my(R,  a);
  R = (2^size+1) * y^(poldegree(pol_field)-1);
  a =  Mod(random(R), pol_field);
  return(a);
};


/* create equation in nf[X] split over nf */
Equation_creation_nf_split(pol_field, degree_eq, size_eq, {type="uniform", norms=0})= {
  my(R, eq, s, b, v, N);
  if (type=="uniform",  R = (2^size_eq+1) * y^(poldegree(pol_field)-1); ,
      s = random(2)*size_eq;
      b = random(2);
      R =  (2^s+1) * y^(poldegree(pol_field)-1);
      );

  v = vector(degree_eq, i, Mod(random(R), pol_field));
  
  if(norms, N = vector(degree_eq, i, norml2(Vec(liftall(v[i])))); );
  eq = 1;
  for (i = 1, degree_eq, eq *= (x - v[i]); );
  if (norms, return([eq, N]); , return(eq) );
};


Equation_creation_nf_multiquad(pol_field, degree_eq, size_roots) = {
  my(v);
  v = vector(degree_eq, i, Random_el_nf(pol_field, size_roots));
  return(prod(i=1, degree_eq, x^2-v[i]));   
};


/* #################### FUNCTIONS SPECIFIC TO RELATIVE EXT #################### */

/* EXTENSION CREATION */

/* generate an irr. pol. over K[X] given
   pol_ground = def. pol of K
   ext_dim = degree of new pol
   size_coeff = size_coeff of new pol */

Pol_ext_creation(pol_ground, ext_dim, size_coeff, {varg=z, vare=y}) =
  {
    my(q, pg);

    pg = subst(pol_ground, variable(pol_ground), varg);
   
    q = Mod(Pol_eq_creation_field(ext_dim-1, size_coeff, pg, vare, varg), pg);
   
    q = q + vare^ext_dim;
   
    while(!polisirreducible(q),

	  q = Mod(Pol_eq_creation_field(ext_dim-1, size_coeff, pg, vare, varg), pg);
   
	  q = q + vare^ext_dim;
   	 
	  );
   
    return(q);
  };


/* create rel ext. fom scratch */
/* dim = vector of dim
   size_coeff = vec of size_coeff  */
Pol_rel_ext_creation(dim, size_coeff, {varg=z, vare=y}) =
  {
    my(p, q);
   
    p = Pol_field_creation(dim[1], size_coeff[1], varg);

    q = Pol_ext_creation(p, dim[2], size_coeff[2], varg, vare);

    return([p, q]);
  };


/* ELEMENTS AND POL CREATION */

/* create g in Ke/Kg ; output is a multivar. pol. */
/* one can retrieve g in vector form by applying Multipol_to_vec */
Rel_element_creation(dim_vec, s_coeff, {varg=z, vare=y})=
  {
    my(a,w);
    w = vecsum(powers(varg, dim_vec[1]-1));
    a = vector(dim_vec[2], i,
	       random((2*2^s_coeff)*varg^(dim_vec[1]-1))-2^s_coeff*w);
    return(Pol(a, vare));
  };


/* create g in Ke/Kg ; output is a multivar. pol. */
/* one can retrieve g in vector form by applying Multipol_to_vec */
Rel_element_creation_mod(pol_vec, s_coeff, {varg=z, vare=y})=
  {
    my(a,w,dg, de);
    dg = poldegree(pol_vec[1]);
    de = poldegree(pol_vec[2]);
   
    w = vecsum(powers(varg, dg-1));
    a = vector(de, i,
	       Mod(random((2*2^s_coeff)*varg^(dg-1))-2^s_coeff*w, pol_vec[1]));
    return(Mod(Pol(a, vare), pol_vec[2]));
  };


/* create Pol in Ke/Kg[x] ; output is a multivar. pol. */
/* one can retrieve Pol in vector form by applying Multipol_to_vec */
Rel_pol_creation(dim_vec, deg, s_coeff, {varp=x, varg=z, vare=y})=
  {
    my(a);
    /* w = vecsum(powers(varp, deg)); */
    a = vector(deg+1, i,
	       Rel_element_creation(dim_vec, s_coeff, varg, vare));
    return(Pol(a, varp));
  };


/* create Pol in Ke/Kg[x] ; output is a multivar. pol. */
/* one can retrieve Pol in vector form by applying Multipol_to_vec */
Rel_pol_creation_mod(pol_vec, deg, s_coeff, {varp=x, varg=z, vare=y})=
  {
    my(a);
    /* w = vecsum(powers(varp, deg)); */
    a = vector(deg+1, i,
	       Rel_element_creation_mod(pol_vec, s_coeff, varg, vare));
    return(Pol(a, varp));
  };



/* create equation of the form h = q(g) with q in ZZ[X] or Ke/Kg[X] */
Rel_equation_creation(pol_vec, s_root, degree_eq, s_eq, type)=
  {
    my(g, h, q, field_dim);
   
    field_dim = poldegree(pol_field);
   
    g = Rel_element_creation_mod(pol_vec, s_root, variable(pol_vec[1]),\
				 variable(pol_vec[2]));
  
    if (type==0,

	q = Pol_creation_integral(degree_eq, s_eq); ,
       
	q = Rel_pol_creation(dim_vec, degree_eq, s_eq, x, variable(pol_vec[1]),\
			     variable(pol_vec[2]));

	);
   
    h = subst(q, x, g);
   
    /* h = Vecrev(lift(h)); */
   
    return([g, h, q])
  };


/* create equation in Ke/Kg[X] split over Ke/Kg */
Rel_equation_creation_split(pol_vec, degree_eq, size_eq)=
  {
    my(v, eq);
   
    v = vector(degree_eq, i, Rel_element_creation_mod(pol_vec, size_eq));
    eq = prod(i = 1, degree_eq, x - v[i]);
   
    return([eq, v])
  };


/* create equation in Ke/Kg[X] split over Ke/Kg */
Rel_equation_creation_multiquad(pol_vec, degree_eq, size_eq)=
  {
    my(v, eq);
   
    v = vector(degree_eq, i, Rel_element_creation_mod(pol_vec, size_eq));
    eq = prod(i = 1, degree_eq, x^2 - v[i]);
   
    return([eq, v])
  };


