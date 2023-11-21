# Code in support of arXiv:2305.17425 [math.NT]
# Copyright 2023, Olivier Bernard, Andrea Lesavourey
# GPL-3.0-only (see LICENSE file)

proof.arithmetic(False);   # Essentially, use pseudo-primality tests everywhere

# ----------------------------------------------------------------------------------
# Some efficient (dedicated to cyclotomic fields) number theoretic functions

# Returns next prime = res mod M greater than "next"
def next_prime_mod_m(M, res, next=2**24):
    _p0 = ZZ(next);
    _p0 = _p0 - ZZ(_p0.mod(M)) + ZZ(res);
    if (_p0 <= ZZ(next)): # Could even test only == 
        _p0 += M;
    while not is_prime(_p0):
        _p0 += M;
    return _p0;


# Primes above p in K
def nf_primes_above(K, p, __idx_start=0):
    r"""
    compute primes above ideal p : output them as two-elements representation
    """
    assert(is_prime(p));
    
    Fp_y   = PolynomialRing(GF(p), name='y');      y = Fp_y.gen();
    Z_x    = PolynomialRing(Integers(), name='x'); x = Z_x.gen();
    f      = Fp_y(K.defining_polynomial());
    f_fact = f.factor();
    return [(p, Z_x(_fact[0])) for _fact in f_fact]


# Does the m-th cyclotomic number field contain inert primes ?
def cf_can_inert(m):
    r"""
    compute whether one can find an inert prime in Q(zeta_m)
    """
    f = factor(m);
    if len(f)>2:
        return False;
    elif len(f)==2 and (f[0][0] > 2 or f[0][1] > 1):
        return False;
    elif f[0][0]==2 and f[0][1] > 2:
        return False;
    else:
        return True;


# K.conductor() is very inefficient, in our case we know K.pol() = some cyclotomic polynomial
def cf_conductor(K):
    return K.gen().multiplicative_order();

# ----------------------------------------------------------------------------------
# Subproduct trees (recursive implementation)

# Output is a nested list of:
#     [ root, left tree, right tree ],
# where root contains the product of the roots of the subtrees.
# A leaf (each element of the input list is a leaf) is [root] and has length 1 instead of 3.
def product_tree(L):
    if len(L) == 1:
        return L;

    mid   = len(L)//2;
    left  = product_tree(L[:mid]);
    right = product_tree(L[mid:]);
    root  = left[0]*right[0];

    return [root, left, right];


# Compute modular values of u modulo the leafs of the given product tree
def mod_tree(T, u):
    mod_u = u.mod(T[0]);
    if len(T) == 1:
        return [ mod_u ];
    else:
        return mod_tree(T[1], mod_u) + mod_tree(T[2], mod_u);


# Same, but for a list of values to be reduced.
# For each leaf of index i, out[i] is the list all_u mod leaf[i].
def mod_tree_all(T, all_u):
    mod_all_u = [ _u.mod(T[0]) for _u in all_u ];
    if len(T) == 1:
        return [ mod_all_u ];
    else:
        return mod_tree_all(T[1], mod_all_u) + mod_tree_all(T[2], mod_all_u);



# ----------------------------------------------------------------------------------
# p-adic lifting (Hensel)

# get suitable prime for p-adic lifting
def cf_padic_get_inert(K, e):
    m = K.conductor()
    n = euler_phi(m)
    
    pp, ee = is_prime_power(m, get_data=True); # Since 1 < m != 2[4], only way is m=4 or pp^ee
    # An element of order phi(m) ?
    r = ZZ(IntegerModRing(m).multiplicative_generator());
    # => inert primes are = r [m] (many choices)

    p = randint(2**60, 2**61)
    p = next_prime_mod_m(m, r, next=p); q = ZZ(p**n);

    if ee==1: # 1/e does exist mod (p^n - 1)/e
        _d = gcd(e,(q-1)//e);
        while _d != 1:
            p = next_prime_mod_m(m, r, next=p); q = ZZ(p**K.degree()); 
            _d = gcd(e,(q-1)//e);
    return p, q


# We compute the "inverse" e-th root: a^(1/e - 1)
# So Let g(x) = a^{e-1} x^e - 1
#        xi+1 = xi - (1/e).xi.g(xi)
# --> Converges to a^{1/e-1}
# --> a.xk = a^{1/e} mod p^k
def cf_padic_lifting(elts, exps, e, B):
    Zx.<x> = PolynomialRing(ZZ);
    K = elts[0].parent();
    mK = cf_conductor(K);
    assert(cf_can_inert(mK)); # Needs inert primes (i.e. cyclic ZZ/mZZ*)

    pp, ee = is_prime_power(m, get_data=True);
    p, q = cf_padic_get_inert(K, e);

    # pp, ee = is_prime_power(m, get_data=True); # Since 1 < m != 2[4], only way is m=4 or pp^ee
    # # An element of order phi(m) ?
    # r = ZZ(IntegerModRing(m).multiplicative_generator());
    # # => inert primes are = r [m] (many choices)
    # _d = 0
    # p = randint(2**60, 2**61)
    # if ee==1: # 1/e does exist mod (p^n - 1)/e
    #     while _d!=1:
    #         p = next_prime_mod_m(m, r, next=p);   print("inert p: ", p, "m: ", m, "r: ", r);
    #         Fp = GF(p, proof=False);
    #         q = ZZ(p**K.degree()); Fq = GF(q, proof=False, modulus=K.polynomial(), name='u'); # if modulus not given, takes forever
    #         _d = gcd(e,(q-1)//e); 
    # else:     # 1/e cannot exist mod (p^n - 1)/e
    #     while _d==0:
    #         p = next_prime_mod_m(m, r, next=p);   print("inert p: ", p, "m: ", m, "r: ", r);
    #         Fp = GF(p, proof=False);
    #         q = ZZ(p**K.degree()); Fq = GF(q, proof=False, modulus=K.polynomial(), name='u'); # if modulus not given, takes forever
    #         _d = gcd(e,(q-1)//e);

    Fp = GF(p, proof=False);
    Fq = GF(q, proof=False, modulus=K.polynomial(), name='u'); # if modulus not given, takes forever
    print("ff done", flush=True);
    Fpy.<y> = PolynomialRing(Fp); gp = Fpy(K.polynomial());
    Fqz.<z> = PolynomialRing(Fq);
    print("poly done",flush=True);

    # Hensel's height
    log2_k = ceil(max(log(ln(2*B)/ln(p), 2),0)); # [Warn] Protection against small cases giving negative log
    k      = 2**log2_k;
    # print("Hensel's height: {} iterations".format(log2_k));
    Rpk     = IntegerModRing(p**k);
    Rpk.<t> = PolynomialRing(Rpk);
    gt      = Rpk(K.polynomial());
    
    # Prod mod p^k
    t = cputime();
    se_Pk = prod( power_mod(Rpk(_u.polynomial()), _e, gt) for _u,_e in zip(elts,exps) ).mod(gt);
    a_Pk  = power_mod(se_Pk, e-1, gt); # Say, a = se^{e-1}
    t = cputime(t);
    # print("time se^{{e-1}} mod P^k: {:.2f}".format(t)); 
    # print("se^{{e-1}} mod P^k:", a_Pk);
    
    # Root in base residue field: first approximation
    t = cputime(); a_P = Fq(Fpy(a_Pk)); t = cputime(t);
    if ee != 1:
        # We should use [AMM77] generalization of Tonnelli-Shanks...
        # Using Cantor-Zassenhaus instead. Seems good enough ?
        pol_pow = a_P*z**e - 1; # [NB] Should compute this mod gt (if e is about not small)
        t = cputime(); is_P  = pol_pow.roots()[0][0]; t = cputime(t); print("roots: {:.2f}".format(t));
    else:
        inv_e = ZZ(-mod(1/e, (q-1)//e)); # careful, this is not mod(-1/e,(q-1)//e) if e not prime
        is_P  = a_P**inv_e;
    # assert(is_P**e * a_P == 1); # Takes time
    # print("a^{{1/e-1}}:", is_P);
    
    # Want to lift root mod P to root mod P^2 etc.
    is_P = is_P.polynomial();
    _k   = 1; # Just in case log2_k = 0
    for _i in range(log2_k):
        # print("\n"+"-"*15+"Lift stage {}".format(_i)+"-"*15);
        _k = 2**(_i+1);
        Fpky.<yk> = PolynomialRing(IntegerModRing(p**_k));
        gp  = Fpky(gt); # Fpky(K.polynomial()); 
        is_P = Fpky(is_P);
        a_P  = Fpky(a_Pk);
        # print("old sol: {}".format(is_P));
        # xi+1 = xi - 1/e xi f(xi)
        eval_isP = a_P * power_mod(is_P, e, gp) - 1;
        # print("f(xi) mod p^{} = {}".format(_k, eval_isP));
        new_isP = is_P - mod(1/e, p**_k) * (is_P * eval_isP).mod(gp);
        # print("xi+1 = ", new_isP);
        # assert( (power_mod(new_isP, e, gp) * a_P).mod(gp) == 1 ); # Takes time
        is_P = new_isP;
    assert(_k == k);
    s_P = (Rpk(is_P)*se_Pk).mod(gt);
    
    print("\n"+"-"*15+"End of p-adic lifting" + "-"*15);
    # assert(power_mod(s_P, e, gt) == se_Pk); # Takes time
    s_P = Zx(s_P);
    s   = K([ (_s if _s < p**k/2. else _s - p**k) for _s in s_P ]);
    # print("SOLUTION: ", s);

    return s;



# ----------------------------------------------------------------------------------
# p-adic reconstruction
# [See] Pohst / Fieker / Belabas / Roblot

# Evaluation of rmax according to [Bel03, Lem.3.8]
def __eval_rmax(La):
    _d = La.nrows();
    GLa,PLa = gram_schmidt_ortho(La);
    iPLa = rc_mat_inverse(PLa);
    rmax = min( 1/2/sqrt(sum(iPLa[_j][_i]**2/(GLa[_j].norm())**2 for _j in range(_d))) for _i in range(_d));
    return rmax;


def belabas(elts, exps, e, lurs=[]):
    red_fct = lll_ZZ;
    K = elts[0].parent();
    la_se = logarg_mpow(lurs, exps);
    la_s  = logarg_pow(la_se, 1/e); # Any argument is fine for pinf[0] (/!\ not consistent with pinf[i])
    m    = cf_conductor(K);
    n    = K.degree();

    if cf_can_inert(m): # call p-adic lifting
        print("can use hensel")
        Be = sum(exps[j]* log(vector(ComplexField(300), elts[j].complex_embeddings(prec=300)).norm(infinity),2) for j in range(len(elts)));
        # print("log h(s^e) = {:.2f}".format(Be));
        B  = Be / e;
        # print("log h(s) = {:.2f}".format(B));
        B =  2**B;
        return cf_padic_lifting(elts, exps, e, B)

    
    # Eval prec
    b_sol  = log(logarg_t2_norm_cf(la_s),2); print("log2 t2(sol): {:.2f}".format(b_sol));
    # Find biggest possible inertia degree
    Zm = IntegerModRing(m);
    fp = max([ _a.multiplicative_order() for _a in Zm if _a.is_unit()]);
    # print("Inertia degree fp={}".format(fp));
    resp = [ _a for _a in Zm if _a.is_unit() and _a.multiplicative_order() == fp ]; assert(len(resp)>0); resp=resp[0];
    p = next_prime_mod_m(m, resp, next=m);
    # print("some prime : {}, fp={}".format(p,fp));
    # ---> Consider Lenstra bound on rmax ([Bel03, end p.17]) and plug in min bi ~ N(P^a)^{1/n}.(1.022)
    #      Condition rmax > min bi / .. > coeff_norm(x) yields:
    #         a > n/ln(N(p)) (ln |x| + ln 2 + (n(n-1)/4-1).ln g_lll ) [ NB: for cyclo, |x| ~ t2(x)/sqrt(n) ]
    g_lll  = 1.022; # Could use a more optimistic value, eg. 1.016
    p_prec = ceil(n/fp/ln(p)*(ln(2)*(1+b_sol)-1/2*ln(n)+(n*(n-1)/4.-1)*ln(g_lll)));
    # next power of 2 is faster to compute (pid^a, lifting).
    print("1st eval: {}".format(p_prec));
    p_prec = 2**(ceil(log(p_prec,2.)));
    print("1st round mod p^a, a={}".format(p_prec));
    # Compute LLL HNF(p^a)
    pid = K.ideal(nf_primes_above(K,p)[0]);
    # [Sage] Costly version in Sage
    # print("pid:", pid);
    # t = cputime(); pid_a = pid**p_prec; t = cputime(t); print("Calcul pid^{}: {:.2f}".format(p_prec, t));
    # t = cputime(); HP_A  = matrix(ZZ, [_a.list() for _a in pid_a.basis()]); t = cputime(t); HP_A  = row_lower_hermite_form(HP_A);
    # print("mat(pid^{}): {:.2f}".format(p_prec, t));
    # //--- [END Sage]
    # [Pari] We want to use Pari for the computation of pid^a, otherwise it takes forever
    Kb = pari.nfinit(K);
    ppid = pari(pid);
    t = cputime(); HP_a = Kb.idealpow(ppid, p_prec); HP_A = matrix(HP_a).transpose(); t = cputime(t);
    print("Calcul pid^{} + HNF: {:.2f}".format(p_prec, t));
    # Convert to HNF in z^i basis
    t = cputime();
    HP_A = row_lower_hermite_form(matrix(ZZ,[ vector(K(Kb.nfbasistoalg(pari(_u).Col()))) for _u in HP_A ]))
    t = cputime(t);
    print("Convert to inc basis {:.2f}".format(t));
    # //--- [END Pari]

    # LLL Reduction of the HNF
    t = walltime(); L_A, _ = red_fct(HP_A); t = walltime(t); print("LLL: {:.2f}".format(t));
    rmax = __eval_rmax(L_A);
    print("ln rmax: {:.2f}, ln target: {:.2f}".format(ln(rmax), ln(RealField()(2**(b_sol)/sqrt(n)))));
    while (rmax < 2**(b_sol)/sqrt(n)):
        p_prec = p_prec*2;
        print("--> increase a={}".format(p_prec));
        # [Sage] version
        # pid_a  = pid_a^2;
        # HP_A = row_lower_hermite_form(matrix(ZZ, [_a.list() for _a in pid_a.basis()]));
        # [Pari] version
        HP_a = Kb.idealpow(HP_a, 2); HP_A = matrix(HP_a).transpose();
        HP_A = row_lower_hermite_form(matrix(ZZ,[ vector(K(Kb.nfbasistoalg(pari(_u).Col()))) for _u in HP_A ]));
        # //--- [END Pari]
        L_A, _ = red_fct(HP_A);
        rmax = __eval_rmax(L_A);
        print("ln rmax: {:.2f}, ln target: {:.2f}".format(ln(rmax), ln(RealField()(2**(b_sol)/sqrt(n)))));
    # Target
    # 1. Obtain s mod P
    # 2. Obtain s mod P^a or (s mod p) mod P^a , by Hensel lifting
    # 3. CVP(P^a, s mod P^a) gives s (or s mod p)
    Fp = GF(p, proof=False); Fpy.<yp> = PolynomialRing(Fp); gp = Fpy((pid.gens()[1]).polynomial());
    q  = p**fp; Fq = GF(q, proof=False, modulus=gp, name='u');
    se_p = prod( Fq(Fpy(_u.polynomial()).mod(gp))**_e for _u,_e in zip(elts,exps) );
    ae_p = se_p**(e-1);
    # Find an e-th root in Fq:
    if (gcd((q-1)/e, e) == 1): # see [BV04] Efficient computation of roots in finite fields
        assert(e.divides(q-1));
        inv_e = ZZ(-mod(1/e,(q-1)/e)); # [Warn] if e is not prime, mod(-1/e, (q-1)/e) != - mod(1/e, (q-1)/e)
        # We want the latter
        x0   = ae_p**inv_e; # x0 = prod^{(1-e)/e} mod Pid
    else: # For now, resort to factorization algo
        # Could be costly, should try to use that we need only one solution to cut the exploration tree in Cantor-Z
        # Or use [AMM77] On taking roots in finite fields
        print("Problem gcd((q-1)/e,e) > 1");
        Fqy.<yq> = PolynomialRing(Fq);
        t=cputime(); x0   = (yq**e * ae_p - 1).roots()[0][0]; t=cputime(t);
        print("Roots mod q using Coppersmith in {:.2f}".format(t));
        # assert(x0**e * ae_p == 1), "Fq:{}\nsp:{}\nx0:{}\nae_p:{}".format(Fq.modulus(),se_p,x0,ae_p);

    # Need product at precision p_prec
    Rpk = IntegerModRing(p**p_prec); Rpk.<t> = PolynomialRing(Rpk);
    gt  = Rpk(HP_A[fp].list()); assert(Fpy(gt) == gp);
    se_pk  = prod( power_mod(Rpk(_u.polynomial()), _e, gt) for _u,_e in zip(elts,exps) ).mod(gt);
    ae_pk  = power_mod(se_pk, e-1, gt); # NB:  Fpy(ae_pk) !=  ae_p, since P^a \not| p

    # A small Hensel lifting to prec a
    x0 = x0.polynomial();
    for _i in range(ceil(log(p_prec,2))):
        _k = min(2**(_i+1), p_prec); print("Lift to prec {}".format(_k));
        Fpky.<yk> = PolynomialRing(IntegerModRing(p**_k));
        gp = Fpky(gt);
        x0 = Fpky(x0);
        ae_p = Fpky(ae_pk);
        eval_f = ae_p*power_mod(x0,e,gp)-1;
        xi = (x0 - mod(1/e,p**_k)*(x0 * eval_f)).mod(gp);
        # assert( (power_mod(xi,e,gp)*ae_p).mod(gp) == 1);
        x0 = xi;
        xpa = (Rpk(x0)*se_pk).mod(gt);
        # print("sol mod P^a ? ", xpa); # So far, seems good

    # [Warn] CVP needs to compute the GSO of L_A at right precision
    v,y=cvp_babai_NP(L_A, vector(ZZ,K(xpa).list()));
    s = K((vector(K(xpa).list()) - v).list());
    return s;



# ----------------------------------------------------------------------------------
# CRT-based method

def root_mod_p(elts, exps, p, eth):
    assert(gcd(p-1, eth)==1); # Otherwise we have multiple choices

    K = elts[0].parent();
    # Define Fp/[y]/1/e mod p
    Fp  = GF(p);
    Fpy = PolynomialRing(Fp, name='y');
    inv_e = ZZ(mod(1/eth, p-1));
    # Orbit above p, plus CRT stuff
    orb_p = nf_primes_above(K, p);
    gx_p  = [ Fpy(_pid[1]) for _pid in orb_p ];
    
    # Subtree of products
    t = cputime(); gx_tree = product_tree(gx_p); t = cputime(t); print("Product tree: {:.2f}".format(t));
    # print(gx_tree); # NB: it contains the full product, we could remove that.
    
    # Residue of the root mod each pid
    t=cputime(); su_p = [ Fpy(_su.polynomial()) for _su in elts ]; t=cputime(t);
    print("Map in Fp[x]: {:.2f}".format(t));
    t=cputime(); su_gx_p = [ [Fp(x) for x in _l] for _l in mod_tree_all(gx_tree, su_p)]; t=cputime(t);
    print("All residues: {:.2f}".format(t));

    # Accumulate and roots.
    res_idp = [];
    T = cputime();
    for p_idx in range(len(gx_p)):
        t = cputime(); _res_se    = prod(su_gx_p[p_idx][j]**ZZ(mod(exps[j],p-1)) for j in range(len(elts))); # print("\tproduct: {:.2f}".format(cputime(t)));
        _res_s     = _res_se**inv_e;
        # assert(_res_s^eth == _res_se);
        res_idp.append(Fpy(_res_s));
        # -- End for
    print("All sols mod P|p: {:.2f}".format(cputime(T)));
    # CRT dans K --> sol mod p
    T = cputime();
    res_p = CRT(res_idp, gx_p);
    print("CRT mod p: {:.2f}".format(cputime(T)));
    return res_p;


# Computation of e-th roots through CRT in number field L
# Q(zeta_e) should not be a subfield of L
def cf_root_crt(elts, exps, e, L, m):
    r"""
    Computes B^(1/e) in OL where B = prod U[i] ^ Exps[i]
    
    Input:
    - 'elts': number field elements; belongs to L
    - 'exps': list of integers; prod elts[i]^exps[i] is of the form A^e
    - 'e': integer; e is the exponent
    - 'L': number field; field where we want to retrieve a e-th root;
    Q(zeta_e) is *not* a subfield of L
    """
    assert(gcd(m,e) == 1); # At the moment, only works with split primes p=1[m], and gcd(p-1,e)==1
    Zx.<x> = PolynomialRing(ZZ);
    bound = sum([abs(ZZ(mod(exps[i],e))) * abs(log(vector(elts[i]).norm(2), 2))
                 for i in range(len(elts))]); # print(bound);
    bound = 2**(1.65*(bound/e).round());
    crt_p = [];
    _p0 = next_prime(2**60);
    _p0 = _p0 + m - _p0.mod(m) + 1; # = 1 mod m
    while (prod(crt_p) < 2 * bound): # factor 2 because of \pm signs
        while ((not is_prime(_p0)) or gcd(_p0-1, e) > 1): # gcd > 1 better than _p0%e==1 when e is a prime-power
            _p0 += m;
        crt_p.append(_p0);
        _p0 += m;
    Pp = prod(crt_p); # we need it to map into -P/2, P/2

    t = cputime(); crt_root = [ Zx(root_mod_p(elts, exps, _p, e)) for _p in crt_p ]; t = cputime(t);
    t = cputime();
    root_big = CRT(crt_root, crt_p) if len(crt_p)>1 else crt_root[0]; # Weird Sage bug if lists have only one element
    root = L([(_s if _s < Pp/2. else _s - Pp) for _s in root_big ]);  # Mapping into -P/2, P/2
    t = cputime(t);

    return root;



# ----------------------------------------------------------------------------------
# Generalized Couveignes' method

# Choice of primes: totally split in K/Q, inert in L/K
# This means p = 1 mod mK, p = (any gen of Z/mLZ*) mod mL/mK ---> congruence condition mod mL
# If possible, we want p != 1 [e^2] to ease the e-th root in OL/p
def cf_couveignes_good_prime(L, K, m, next=1, e=1):
    r"""
    Find a good prime for couveignes' method relative to L / K

    INPUT:
    - ''L'' -- cyclotomic number field;
    - ''K'' -- cyclotomic number field; K is a subfield of L
    - ''m'' -- [mK, mL] conductors of K and L respectively
    - ''next'' -- return a prime (strictly) greater than next
    - ''e'' -- this is for computing an e-th root (prefer p != 1 [e^2])

    OUTPUT: a tuple p, Ip, Jp such that
    - p is a prime interger such that all prime ideals of K above it are inert
      in L/K
    - Ip is the list of prime ideals of K above p
    - Jp is the list of prime ideals of L above p
    """
    mK = m[0];
    mL = m[1];
    mLK = mL / mK;
    
    # Choose primes p = (a mod mL/mK, 1 mod mK) mod mL
    a   = ZZ(IntegerModRing(mLK).multiplicative_generator());
    res = CRT([1,a],[mK,mLK]); # target residue mod mL: split in Q(zmK), inert mod mLK)
    p = next_prime_mod_m(mL, res, next=next);
    if (e != 1) and (mod(mK, e**2) != 0):
        while (mod(p-1, e**2) == 0): # For e-th roots of relative norms, it might be more comfortable if p!=1[e^2]
            p = next_prime_mod_m(mL,res,next=p);
    
    Ip = nf_primes_above(K, p);
    Jp = nf_primes_above(L, p);
    len_K = len(Ip); assert(len_K == K.degree());
    len_L = len(Jp); assert(len_L == len_K);

    Fpy.<y> = PolynomialRing(GF(p));
    Jp1 = [];
    for _ip in Ip:
        _gK = _ip[1];
        _gLK = _gK(L.gen()**mLK).polynomial(); # generators of ideals in L
        for i in range(len(Jp)):
            _jp = Jp[i];
            _gL = _jp[1];
            if Fpy(_gLK.mod(_gL)) == 0:
                Jp1 += [_jp];
                _a = Jp.pop(i);
                break;
    Jp = Jp1;
    return p, Ip, Jp;


# Relative norm maps in extension of cyclotomic fields
def cf_relative_norms(A, L, K):
    r"""
    Computes N_{L/K}(a) for a in A
    
    Input:
    - 'A': number field elements; belong to L
    - 'L': number field;
    - 'K': number field; subfield of L verifying good properties
    ( gcd([L:K], e) == 1, zeta_e belongs to K)
    """

    Zx.<x> = PolynomialRing(ZZ);
    Ky.<y> = PolynomialRing(K);
    Lz.<z> = PolynomialRing(L);

    mK = cf_conductor(K)
    mL = cf_conductor(L)
    mLK = ZZ(mL/mK); # Raises an error if by accident mK \not| mL

    # Construct embeddings K --> L: sigma_s with s mod mLK invertible, s mod 3mK= 1
    embLK = [ CRT([i,1],[mLK,mK]) for i in range(1,mLK) if gcd(i,mLK)==1 ];
    ext_pol = prod(z - L.gen()**i for i in embLK);
    ext_pol = Ky( { _key: K(_val) for (_key,_val) in ext_pol.dict().items()} );
    # print("ext_pol", ext_pol);
    
    # Computation of relative norms by computing the product of conjugates
    xml_pol = x**mL-1;
    phiL = Zx(L.defining_polynomial());
    # [Slow] Compact version, but we need to apply the modulo at each step of the product
    # NA = [ K(list(prod(Zx(_a.polynomial())(x**i).mod(x**mL-1) for i in embLK))[::mLK]) for _a in A];
    NA = [ ];
    t = cputime();
    for _a in A:
        _na = Zx(1);
        _ax = _a.list();
        in_t = cputime();
        for _i in embLK: # Necessary to apply mod x^mL-1 at each step
            _axi = Zx(sum([_ax[_k]*x**ZZ(mod(_k*_i,mL)) for _k in range(len(_ax))]));
            _na  = (_na*_axi).mod(xml_pol);

        in_t = cputime(in_t); # print("inner ", in_t);
        _na = _na.mod(phiL);
        _na = _na(y).mod(ext_pol); 
        assert(_na.degree() == 0);
        
        NA.append(K(_na));
    t = cputime(t);
    print("all relative norms: {:.2f} [{:.2f}/elt]".format(t, t/len(A)));

    # [Test] Resultants --------------------------------------------------------------
    # t = cputime();
    # ZZut.<u,v> = PolynomialRing(ZZ,2);
    # Ext_pol = [ _c.polynomial()(v) for _c in ext_pol.list() ]; assert(len(Ext_pol) == euler_phi(mLK)+1);
    # Ext_pol = ZZut(sum( u**_i*Ext_pol[_i] for _i in range(len(Ext_pol)) ));
    # NA2 = [ ZZut(_a.polynomial()(u)).resultant(Ext_pol,u).mod(K.defining_polynomial()(v)) for _a in A ];
    # NA2 = [ K(_na(1,x)) for _na in NA2 ];
    # t = cputime(t);
    # print("all relative norms using resultants: {:.2f} [{:.2f}/elt]".format(t, t/len(A)));
    # assert(NA == NA2);
    # //--- [END Resultants] Only works for small values, and always twice slower anyway. Does not scale at all.

    return NA;


# Local method using Couveignes' criterion for one prime **in cyclotomic fields**
def cf_couveignes_mod_p(U, Exps, e, NA, L, K, m, p, PI_p):
    r"""
    Computes B^(1/e) in OL/(p), where B = prod_i U[i]^Exps[i]

    Input:
    - 'U', 'Exps': number field elements and exponents; belong to L and product is of the form A^e
    - 'e': integer; e is the order of the seeked root
    - 'NA': number field element; norm of A relative to L/K times some power of zeta_e,
            where zeta_e is a primitive e-th root of unity in K
    - 'L': number field; field where we want to retrieve a e-th root of B
    - 'K': number field; subfield of L verifying good properties (gcd([L:K], e) == 1, zeta_e belongs to K)
    - 'm': [mK, mL] respective conductors of K, L
    - 'p', integer; computations done modulo p
    - 'PI_p': [Ip,Jp] Ip are prime ideals of K above p, Jp are prime ideals of L above p (same order)
    """
    mK = m[0];
    mL = m[1];
    mLK = ZZ(mL / mK);
    
    Fp  = GF(p, proof=False);
    Fpy = PolynomialRing(Fp, name='y'); y = Fpy.gen();

    gx_p_K = [Fpy(_ip[1]) for _ip in PI_p[0]]  # generators of ideals in K
    gx_p_L = [Fpy(_jp[1]) for _jp in PI_p[1]]  # generators of ideals in L
    
    # residue fields
    q  = p**gx_p_L[0].degree();
    Lp = [GF(q, modulus=_gx, name='A', proof=False) for _gx in gx_p_L];

    # embed norm in residue fields
    NAp = [Fp (Fpy(NA.polynomial()).mod(_gx)) for _gx in gx_p_K];

    gt = Fpy(L.polynomial())
    BL = prod( power_mod(   Fpy(_u.polynomial()), _e, gt ) for _u,_e in zip(U,Exps) ).mod(gt);
    BLp = [ _lp(BL) for _lp in Lp ];

    my_t = cputime()
    assert(ZZ(e).divides(p-1));
    Fp_gen = Fp.multiplicative_generator();
    ze_p   = [ Fp_gen**(_k*(p-1)/e) for _k in range(e) ];
    # print("time for ze mod p: ", cputime(my_t));

    # e-th roots of embeddings of B in residue fields
    my_t = cputime();
    if (gcd(e, (p-1)//e) == 1):
        print("e-th root in Fq using inversions");
        inv_e  = ZZ(mod(1/e, (q-1)/e));
        BLp_e  = [ _blp**inv_e for _blp in BLp ]; # Could just apply ^1/e ***before*** mapping in Lp
    else:
        print("WARN: e-th root in Fq using generic methods", end="\t");
        BLp_e = [];
        if e < 10:         # factorisation of polynomials z^e - A in residue fields if e is small
            print("Calling eq.roots()");
            for i in range(len(BLp)):
                Q.<z> = PolynomialRing(Lp[i])
                eq = z**e - BLp[i]
                pr = eq.roots()
                pr = pr[0][0]#[_pr[0] for _pr in pr]
                BLp_e += [pr]
        else:
            print("Calling nth_root");
            BLp_e = [_blp.nth_root(e, all=False) for _blp in BLp];
    # print("time taken for eth roots is: ", cputime(my_t));

    # norms of previous roots
    my_t = cputime();
    frob  = (q-1)//(p-1); 
    NB_mod_ze = [ Fp(_blp_e**frob)/_na for _blp_e, _na in zip(BLp_e, NAp) ];
    # print("time for relative norms mod ze in Fp: ", cputime(my_t));
    
    # for each residue field, test which norm is good
    ze_pow_idx = [ ze_p.index(_ze) for _ze in NB_mod_ze ];
    # print("time afer index: ", cputime(my_t));
    inv_nLK = mod(1/euler_phi(mLK),e); # normally guaranteed to be invertible
    ALp = [ _blp*ze_p[( -_idx*inv_nLK )] for _blp, _idx in zip(BLp_e,ze_pow_idx) ];

    # reconstruction via CRT
    Ap = CRT([Fpy(_alp.polynomial()) for _alp in ALp], gx_p_L);
    # print("time after CRT: ", cputime(my_t));
    
    return Ap;



# -------------------------------------------------------------------------------------
# The main recursive function
def cf_roots_rec(U, Exps, e, L, mL, size_primes=60):
    r"""
    Computes B^(1/e) in OL where B = prod U[i] ^ Exps[i]
    Recursive version
    
    Input:
    - 'U': number field elements; belongs to L
    - 'Exps': list of integers; prod U[i]^Exps[i] is of the form A^e
    - 'e': integer; e is the exponent
    - 'L': number field; field where we want to retrieve an e-th root of B
    """

    if gcd(mL,e) == 1:               # good case : we can do a CRT à la Thomé
        print("can do crt");
        return cf_root_crt(U, Exps, e, L, mL), 0;
    else:                       # bad case : factor the conductor to know what to do
        if cf_can_inert(mL): # len(fL)==1:          # mL=e^k for some k : we can use plain Hensel
            print("can use hensel");
            # Fits well with real measurements for cyclotomics
            Be = sum(Exps[j]* log(vector(ComplexField(300), U[j].complex_embeddings(prec=300)).norm(infinity),2)
                     for j in range(len(U)));
            # print("log h(s^e) = {:.2f}".format(Be));
            B  = Be / e;
            # print("log h(s) = {:.2f}".format(B));
            B =  2**B;
            return cf_padic_lifting(U, Exps, e, B), 0;
        else:                   # can we use couveignes then ?
            fK = [];
            mK = 1;
            fL = factor(mL);
            for _f in fL:
                if (_f[0]==e) or (euler_phi(_f[0])%e == 0):
                    fK.append(_f);
                    mK *= _f[0]**_f[1];
            mLK = mL // mK;
            # print(mLK)
            if mK==mL or (not cf_can_inert(mLK)):
                # /!\ WORST CASE : we use p-adic reconstruction by [Bel03]
                print("using p-adic reconstruction")
                phi = get_inf_places(L,300);
                lurs = [ logarg_set(_u, phi) for _u in U ];
                return belabas(U, Exps, e, lurs), 0;
            else: # start Couveignes' process
                print("using couveignes' method");
                K = CyclotomicField(mK);
                m = [mK,mL];
                Zx.<x> = PolynomialRing(ZZ);
                nL = L.degree();
                nK = K.degree();
                nLK = ZZ(nL / nK);
                
                # Compute bounds
                my_t = cputime();
                Re = RealField(300);
                bound_u = [abs(log(vector(Re,_u).norm(2), 2)) for _u in U];
                bound_basis = max(bound_u); # much faster if not in symb ring
                # bound_norm = round(2**(bound_basis*nLK)/2 +  log(sqrt(Re(nL))*nLK, 2) /2)
                print("calcul des bornes: {:.2f}".format(cputime(my_t)));
                
                # Compute all relative norms
                my_t = cputime();
                NU = cf_relative_norms(U,L,K);
                print("relative norm computed in: ", cputime(my_t));
                t_norm = cputime(my_t);

                # e-th Root of relative norm of the product
                my_t = cputime();
                NA, _ = cf_roots_rec(NU, Exps, e, K, mK); # Recursive e-th root of relative norm
                print("root of norm computed in: ", cputime(my_t));
                
                # now will start CRT process
                A_crt = [];
                my_t = cputime();
                bound = sum([ abs(_e)*_bu for _e, _bu in zip(Exps,bound_u) ]);
                bound = 2**(1.3*(bound/e).round());
                print("bounds for couveignes: ", cputime(my_t));
                p = randint(2**size_primes, 2**(size_primes+1));
                # p = p- p%mK + 1
                crt_Ip = [];
                crt_p = [];
                crt_Jp = [];
                crt_Ap = [];
                prod_crt_p = 1;
                
                my_t = cputime();
                while prod_crt_p < bound: # computation of suitable primes
                    p, Ip, Jp = cf_couveignes_good_prime(L, K, m, next=p, e=e);
                    prod_crt_p *= p;
                    crt_p.append(p);
                    crt_Ip.append(Ip);
                    crt_Jp.append(Jp);
                print("primes for couv computed in: ", cputime(my_t));
                for i in range(len(crt_p)): # now we'll find the root mod p for each p
                    my_int = cputime();
                    print("treating ", crt_p[i], flush=True);
                    crt_Ap += [cf_couveignes_mod_p(U, Exps, e, NA, L, K, m, crt_p[i],
                                                   [crt_Ip[i], crt_Jp[i]])];
                    print("one couveignes mod p computed in: ", cputime(my_int));
                crt_Ap = [Zx(_ap) for _ap in crt_Ap];

                my_t = cputime();
                A_big = CRT(crt_Ap, crt_p) if len(crt_p)>1 else crt_Ap[0]; # reconstruction through CRT
                A = L([(_a if _a < prod_crt_p/2. else _a - prod_crt_p) for _a in A_big]);
                print("CRT for couveignes in: ", cputime(my_t));
                
                return A, t_norm;


# //-- END OF FILE


# -------------------------------------------------------------------------------------------
# [OBSOLETE] Historic code to erase
# embedding of a cyclotomic field into another
def __unused_cf_elt_embed(g, K, L):
    r"""
    change representation g in K --> g in L, where K and L are both cyclotomic
    fields
    """
    _gx = g.polynomial()
    mK = K.conductor()
    mL = L.conductor()
    mLK = mL / mK
    return L(_gx(L.gen()**mLK))

    
### norm map in finite fields --->  not used anymore
def my_ff_norm(x, q, Q):
    r"""
    compute the norm of x relative to FQ / Fq by applying frob (=generator
    of the relative galois group Gal(FQ / Fq) ) several times
    """
    _e = (Q-1)//(q-1)
    return x**_e

def __unused_my_ff_norm_ext(x, frob, o):
    r"""
    compute the norm of x relative to FQ / Fq by applying frob (=generator
    of the relative galois group Gal(FQ / Fq) ) several times
    """
    q, r = divmod(o, 2);
       
    frob2 = frob**2;
    y = x;
    
    for i in range(q-1):
        y = y**frob2;
        y *= x
    
    tt = cputime()
    y *= y**frob;
    # print("last frob in :", cputime(tt))
    if r == 1:
        y = x * y**frob;

    return y;

def __unused_my_frob(FQ, Fq):
    r"""
    return frob = frobenius relative to FQ / Fq i.e. frob : x -> x^q  ;
    q = p^d and Q = q^n for some d and n
    """
    p = Fq.characteristic();    
    d = log(Fq.order(), p);
    return FQ.frobenius_endomorphism(d);


# [Some notes on Newton methods]
#     - Brent's homepage: https://maths-people.anu.edu.au/~brent/pub/pub226.html
#     - Planetmath:       https://planetmath.org/nthrootbynewtonsmethod
#     - Halley's method, Pade approximant + some insight: https://math.stackexchange.com/questions/1825062/an-apparently-new-method-to-compute-the-nth-root-of-any-complex-number
# 
# 1. Root of f(x) = x^e - a using Newton.
#     f(xi) = 0 mod p^k at precision k
#     xi+1 = xi - f(xi)/f'(xi)
# 2. Suggest computing inverse root first
#      a^{-1/2} --> a^{1/2} = a*a^{-1/2}
#      f(x) = a.x^2 - 1
#    --> No inversion ??
#        "xi+1 = xi - 1/2.xi(ei - 3/4.ei^2)"
#    --> 3rd order convergence
#    --> ei: n, ei^2 : n/3, xi(ei-3/4ei^2): 2n/3
# 3. What does it give for e-th root ?
#        Halley: h(x) = x^((e+1)/2) - a.x^(-(e-1)/2)
#        1/e-1 : g(x) = a*x^(1+1/(e-1)) - 1 --> Pb!! what does mean "^{1+1/(e-1)}" mod p**k ?
#            use g(x) = a^{e-1} x^e - 1
#        Plain : f(x) = x^e - a
# 3.a. Pour Halley: Fe(z) = frac{(e-1) + (e+1)z}{(e+1)+(e-1)z}
#                   xi+1 = xi*Fe(a/xi^e) ---> Need inverse
# 3.b. xi+1 = xi - ((e-1)/e) xi f(xi), puis xn*a = a^{1/e}. 
# 3.c  xi+1 = xi - 1/(ea) xi f(xi),      ---> Need a, inverse of a
#


# Bound for Belabas p-adic reconstruction
# Use as suggested [Bel03, Lem.3.12] without the (3.sqrt(gamma)/2)^{d-1} 
# eval_a = ceil((n*(1/2.*(b_sol-log(n,2))+1))/(fp*log(p,2)));

# Snippet code for p-adic reconstruction if staying in Pari's basis for OK
# print("HNF qu'on veut:", row_lower_hermite_form(matrix(ZZ, [pari.nfinit(K).nfalgtobasis(pari(_a)).list() for _a in pid.basis()])));
# print("HNF Pari", matrix(pari.nfinit(K).idealhnf(ppid)));
# t = cputime(); HP_a = pari.nfinit(K).idealpow(ppid, p_prec); t = cputime(t); print("Calcul pid^{} + HNF: {:.2f}".format(p_prec, t));
# AND
# v,y=cvp_babai_NP(L_A, vector(ZZ, pari(K).nfalgtobasis(xpa).list()));
# s = K (pari(K).nfbasistoalg((vector(pari(K).nfalgtobasis(xpa).list()) - v).list()));
# print("SOLUTION: ", s);


# ---------------------------------------------------------------------
# Relative norms computations with local / global principle ?
# As of now, this is not competitive
def cf_norm_good_primes(L, K, m, q=1):
    r"""
    Find a good prime for norm computation relatively to L / K

    INPUT:

    - ''L'' -- cyclotomic number field;
    
    - ''K'' -- cyclotomic number field; K is a subfield of L
    
    - ''m'' -- integer; m is a product of forbidden primes ; typically the
       conductor of L.

    OUTPUT: a tuple p, Ip, Jp such that

    - p is a prime interger completely split in L
    
    - Ip is the list of prime ideals of K above p
    
    - Jp is the list of prime ideals of L above p
    """
    mK = m[0]
    mL = m[1]
    pL = L.defining_polynomial()
    pK = K.defining_polynomial()
    nL = L.degree()
    mLK = mL / mK

    p = q + mL;

    while (not p.is_prime()) or (mL % p == 0):
        p = p + mL

    Ip = nf_primes_above(K, p);
    Jp = nf_primes_above(L, p);
    len_K = len(Ip)
    len_L = len(Jp)

    my_t = cputime()
    while len_L != nL:
        p = p + mL;
        while (mL % p == 0) or (not p.is_prime()):
            p = p + mL;
        Ip = nf_primes_above(K, p);
        Jp = nf_primes_above(L, p);
        len_K = len(Ip)
        len_L = len(Jp)
    # print("time taken for finding one prime is: ", cputime(my_t))
    
    my_t = cputime()
    Fpy.<y> = PolynomialRing(GF(p))
    Jp1 = []
    for _ip in Ip:
        Jp2 = []
        _gK = _ip[1]
        _gLK = _gK(L.gen()**mLK).polynomial() # generators of ideals in L
        for i in range(len(Jp)):
            _jp = Jp[i]
            _gL = _jp[1]
            if Fpy(_gLK.mod(_gL)) == 0:
                Jp2 += [_jp]
                # _a = Jp.pop(i)
        Jp1 += [Jp2]
        
    
    Jp = Jp1
    return p, Ip, Jp



### relative norm mod one prime and *cyclotomic fields*
def cf_relative_norm_mod_p(a, L, K, m, p, PI_p, zm=0, ze = 0):
    r"""
    Computes relative norm N_{L/K}(a) mod p, where p is given as
    cf_couveignes_good_prime.

    To use when the relative degree [L:K] is large.
       
    
    Input:
    
    - 'a': number field element; belongs to L 
        
    - 'L': number field;
    
    - 'K': number field; subfield of L verifying good properties (gcd([L:K], e) == 1,
      zeta_e belongs to K)
    
    - 'p', integer; computations done modulo p
    
    - 'Ip', list of number field ideals; prime ideals of K above p
    
    - 'Jp', list of number field ideals; prime ideals of L above p
    """
    MY_T = cputime()
    mK = m[0]
    mL = m[1]
    mLK = mL / mK
        
    Fp  = GF(p);
    Fpy = PolynomialRing(Fp, name='y');
    Zx = PolynomialRing(ZZ, name='x');
    if zm == 0:
        zm = L.gen()
    if ze == 0:
        ze = L.gen()**(mL/e)

    Ip = PI_p[0]
    Jp = PI_p[1]

    gx_p_K = [_ip[1] for _ip in Ip]
    gx_p_L = [_jp[1] for _jp in Jp]

    gx_p_K = [Fpy(_gx) for _gx in gx_p_K] # generators of ideals in K
    gx_p_L = [Fpy(_gx) for _gx in gx_p_L] # generators of ideals in L
    
    Kp = [GF(p**(_gx.degree()), modulus=_gx, name='a') for _gx in gx_p_K];
    Lp = [GF(p**(_gx.degree()), modulus=_gx, name='A') for _gx in gx_p_L];

    crt_Na_p = []

    crt_a_p = [Fpy(a.polynomial())(_lp.gen()) for _lp in Lp]; # embed norm in residue fields
    crt_Na_p = [ my_ff_norm(crt_a_p[i],  Kp[i].cardinality(), Lp[i].cardinality())
                 for i in range(len(crt_a_p))]
    Nap = CRT([Fpy(_nap.polynomial()) for _nap in crt_Na_p], gx_p_K);
    
    return Nap;

# OLD VERSIONS ---- TRY NEW VERSION BELOW
# ### relative norm mod one prime and *cyclotomic fields*
# def cf_relative_norms_mod_p(A, L, K, m, p, PI_p, zm=0, ze = 0):
#     r"""
#     Computes relative norm N_{L/K}(a) mod p, where a are elements of list A
#     and p is given as in cf_couveignes_good_prime.

#     To use when the relative degree [L:K] is large.
    
#     Input:
    
#     - 'A': number field elements; belongs to L 
        
#     - 'L': number field;
    
#     - 'K': number field; subfield of L verifying good properties (gcd([L:K], e) == 1,
#       zeta_e belongs to K)
    
#     - 'p', integer; computations done modulo p
    
#     - 'Ip', list of number field ideals; prime ideals of K above p
    
#     - 'Jp', list of number field ideals; prime ideals of L above p
#     """

#     MY_T = cputime()
#     mK = m[0]
#     mL = m[1]
#     mLK = mL / mK
        
#     Fp  = GF(p);
#     Fpy = PolynomialRing(Fp, name='y');
#     Zx = PolynomialRing(ZZ, name='x');
#     if zm == 0:
#         zm = L.gen()
#     if ze == 0:
#         ze = L.gen()**(mL/e)

#     Ip = PI_p[0]
#     Jp = PI_p[1]
#     # JJp = PI_p[1]
    
#     gx_p_K = [_ip[1] for _ip in Ip]
#     gx_p_L = [_jp[1] for _jp in Jp]

#     gx_p_K = [Fpy(_gx) for _gx in gx_p_K] # generators of ideals in K
#     gx_p_L = [Fpy(_gx) for _gx in gx_p_L] # generators of ideals in L
    
#     Kp = [GF(p**(_gx.degree()), modulus=_gx, name='a') for _gx in gx_p_K];
#     Lp = [GF(p**(_gx.degree()), modulus=_gx, name='A') for _gx in gx_p_L];

#     NAp = []

#     for a in A:
#         crt_Na_p = []
        
#         crt_a_p = [Fpy(a.polynomial())(_lp.gen()) for _lp in Lp]; # embed norm in residue fields
#         crt_Na_p = [ my_ff_norm(crt_a_p[i],  Kp[i].cardinality(), Lp[i].cardinality())
#                      for i in range(len(Lp))]
        
#         Nap = CRT([Fpy(_nap.polynomial()) for _nap in crt_Na_p], gx_p_K);

#         NAp += [Nap]

    
#     # # #############   with completely split primes  ###########
#     # crt_Na_p = [[] for _a in A]

#     # for i in range(len(Ip)):
#     #     Jp = JJp[i];
#     #     gx_p_L = [_jp[1] for _jp in Jp]
#     #     gx_p_L = [Fpy(_gx) for _gx in gx_p_L] # generators of ideals in L
        
#     #     my_t = cputime()
#     #     Lp = [GF(p**(_gx.degree()), modulus=_gx, name='A') for _gx in gx_p_L];
#     #     # print("time taken for residue fields is: ", cputime(my_t));


#     #     for j in range(len(A)):
#     #         a=A[j]
#     #         # for a in A:
#     #         crt_a_p = [Fpy(a.polynomial())(_lp.gen()) for _lp in Lp]; # embed norm in residue fields
#     #         crt_Na_p[j] += [prod(crt_a_p)]

#     # for i in range(len(A)):
#     #     Nap = CRT([Fpy(_nap.polynomial()) for _nap in crt_Na_p[i]], gx_p_K);
#     #     NAp += [Nap]
#     #     # NAp += [crt_NA_p[i]]
        
#     return NAp;


### relative norm map in number field
def cf_relative_norms_localglobal(A, L, K, crt_p, crt_Ip, crt_Jp, bound):
    r"""
    Computes N_{L/K}(a)
    
    Input:
    - 'U': number field element; belongs to L
    - 'L': number field;
    - 'K': number field; subfield of L verifying good properties
    ( gcd([L:K], e) == 1, zeta_e belongs to K)
    """

    Zx.<x> = PolynomialRing(ZZ);
    Qy.<y> = PolynomialRing(ZZ);
    my_t = cputime()
    mK = K.conductor()
    mL = L.conductor()
    nL = L.degree()
    m = [mK, mL]

    # print("time taken for technical stuff is: ", cputime(my_t))
    
    nK = K.degree()

    my_t = cputime()    
    
    if crt_p == []:
        p = randint(2**60, 2**61);
        p = p-p%mL + 1
        prod_crt_p = 1;
    else:
        p = crt_p[len(crt_p)-1]
        prod_crt_p = prod(crt_p);
    while (prod_crt_p < bound): 
        p, Ip, Jp = cf_couveignes_good_prime(L, K, m, q=p); # inert seems more efficient
        # p, Ip, Jp = cf_norm_good_primes(L, K, m, q=p);
        prod_crt_p *= p;
        crt_p += [p];
        crt_Ip += [Ip];
        crt_Jp += [Jp];
    
    # print("time taken for computation of primes is:", cputime(my_t));

    NA = []
    NA_crt = [[] for _a in A]
    for i in range(len(crt_p)):
        NAp = cf_relative_norms_mod_p(A, L, K, [mK, mL], crt_p[i], [crt_Ip[i], crt_Jp[i]])
        for j in range(len(A)):
            NA_crt[j] += [NAp[j]]
        # print("time taken for computation of one norm mod p is:", cputime(my_t));

    for i in range(len(A)):
        Na_crt = NA_crt[i];
        Na_crt = [Zx(_nap) for _nap in Na_crt];
        my_t = cputime();
        Na_big = CRT(Na_crt, crt_p) if len(crt_p)>1 else Na_crt[0];
        # print("time taken for CRT is:", cputime(my_t));
        Na = K([(_a if _a < prod_crt_p/2. else _a - prod_crt_p) for _a in Na_big]);
        NA += [Na]

    return NA;



    # gx_p_K = [(pid_fast_gens_two(_pid)[1]) for _pid in Ip] # generators of ideals in K
    # gx_p_L = [(pid_fast_gens_two(_pid)[1]) for _pid in Jp] # generators of ideals in L
    # gx_p_L = [_pid.gens_two()[1].polynomial() for _pid in Jp] # generators of ideals in K
    # gx_p_L = [(_gx(zm**mLK)) for _gx in gx_p_K] # generators of ideals in L
    # gx_p_L = [_gx.polynomial() for _gx in gx_p_L] # generators of ideals in L
    # print("time taken for residue fields is: ", cputime(my_t));
   
    # # embeddings of Kp[i] into Lp[i]
    # Kp = [GF(p**(_gx.degree()), modulus=_gx, name='a') for _gx in gx_p_K];
    # Embp = [FiniteFieldHomomorphism_generic(Hom(Kp[i], Lp[i])) for i in range(len(Lp))];
    # print("after embeddings in:", cputime(MY_T))
    # Zep = [Fpy(ze.polynomial())(_lp.gen()) for _lp in Lp]; # embed ze into residue fieds
    # NAp = [Fpy(NA.polynomial())(_kp.gen()) for _kp in Kp]; # embed norm in residue fields
    # NAp = [Embp[i](NAp[i]) for i in range(len(NAp))]
    # print("after use of emb morph in:", cputime(MY_T))
    # print("after embedding of norm in:", cputime(MY_T))
    
    # Zx = PolynomialRing(ZZ, name='x');
    # BL = prod([L(Zx(Fpy(U[i].polynomial())**Exps[i])) for i in range(len(U))])
    # BLp1 = [_lp(BL.polynomial()) for _lp in Lp]

    # assert(BLp == BLp1);
    # BLp = [ prod( _Lp( Fpy(_u.polynomial()).mod(_gx) )**_e for _u,_e in zip(U,Exps) ) for _Lp,_gx in zip(Lp,gx_p_L) ]
    # BLp = [prod([_lp(Fpy(U[i].polynomial()))**Exps[i] for i in range(len(U))]) for _lp in Lp]
    # print("before computing e-th roots in residue field in:", cputime(MY_T))

    # frob = (q-1)//(p-1);    
    # # norms of previous roots
    # NBp_r = [ [BLp_r[i][j]^frob # Vaguely faster my_ff_norm(BLp_r[i][j], p, q)
    #            for j in range(len(BLp_r[i]))] for i in range(len(Lp))];
    # print("time taken for norms is: ", cputime(my_t));

    # for each residue field, test which norm is good
    # ALp = [];
    # for i in range(len(Lp)):
    #     for j in range(len(NBp_r[i])):
    #         if NBp_r[i][j] == NAp[i]:
    #             ALp += [BLp_r[i][j]]
    #             break;

### Une tentative de faire une fonction couveignes, mais il faut appeler cf_roots qui elle-même appelle Couveignes
### Il faut réfléchir à une archi possible.
### Attention ce code est obsolète il faudra repartir du contenu de cf_roots_rec
### couveignes generalisation in *cyclotomic fields*
def cf_couveignes(U, Exps, e, L, K, m):
    r"""
    Computes B^(1/e) in OL where B = prod U[i] ^ Exps[i]
    
    Input:
    - 'U': number field elements; belongs to L
    - 'Exps': list of integers; prod U[i]^Exps[i] is of the form A^e
    - 'e': integer; e is the exponent
    - 'L': number field; field where we want to retrieve a e-th root of B
    - 'K': number field; subfield of L verifying good properties
    ( gcd([L:K], e) == 1, zeta_e belongs to K)
    """
    Zx.<x> = PolynomialRing(ZZ);
    # Qy.<y> = PolynomialRing(ZZ);
    my_t = cputime()

    #print("time taken for technical stuff is: ", cputime(my_t))
    nL = L.degree()
    zL = L.gen()

    ze = zL**(m[1]/e)
    
    # mK = K.conductor();
    nK = K.degree()
    nLK = ZZ(nL / nK)

    my_t = cputime()

    bound_basis = max([abs(log(vector(_u).norm(2), 2)) for _u in U])
    bound_norm = round(2**(bound_basis*nLK)/2 +  log(sqrt(nL)*nLK, 2) /2)
    
    if nLK < 10:
        NA = prod([(U[i].norm(K))**Exps[i] for i in range(len(U))])
    else:
        NU = cf_relative_norms(U, L, K, [], [], [], bound_norm)
        NA = prod([NU[i]**Exps[i] for i in range(len(U))])

    #print ("time taken for relative norm is: ", cputime(my_t))
    NA = NA.nth_root(e);
    #print ("time taken for root in subfield is: ", cputime(my_t))

    A_crt = [];
    bound = 2**(1.5*(log(max([abs(_b) for _b in B]), 2)/e).round());
    p = 1;
    crt_p = [];
    crt_Ip = [];
    crt_Jp = [];
    crt_Ap = [];
    prod_crt_p = 1;

    my_t = cputime();
    while prod_crt_p < bound:
        p, Ip, Jp = cf_couveignes_good_prime(L, K, m, q=p);
        prod_crt_p *= p;
        crt_p += [p];
        crt_Ip += [Ip];
        crt_Jp += [Jp];
    #print ("time taken for computation of primes is:", cputime(my_t));

    print ([len(_ip) for _ip in crt_Ip])
    for i in range(len(crt_p)):
        # print("prime is: ", crt_p[i]);
        my_t = cputime();
        crt_Ap += [cf_couveignes_mod_p(U, Exps, e, NA, L, K, m, crt_p[i],
                                       [crt_Ip[i], crt_Jp[i]], zL, ze)];
        # print ("time taken for couveignes_mod_p is:", cputime(my_t));
        # print("******************************");
        # print("everything is done")

    crt_Ap = [Zx(_ap) for _ap in crt_Ap];
    my_t = cputime();
    A_big = CRT(crt_Ap, crt_p);
    #print("time taken for CRT is:", cputime(my_t));
    A = L([(_a if _a < prod_crt_p/2. else _a - prod_crt_p) for _a in A_big]);
    return A;


# Relative norm maps in extension of cyclotomic fields
def cf_relative_norms_mod_p(A, L, K, p, k):
    r"""
    Computes N_{L/K}(a) mod p for a in A
    
    Input:
    - 'A': number field elements; belong to L
    - 'L': number field;
    - 'K': number field; subfield of L verifying good properties
    ( gcd([L:K], e) == 1, zeta_e belongs to K)
    - 'p': prime integer
    """
    Zx.<x> = PolynomialRing(ZZ); 
    ZpX.<X> = PolynomialRing(IntegerModRing(p**k));
    Ky.<y> = PolynomialRing(K);
    Lz.<z> = PolynomialRing(L);
    
    time mK = cf_conductor(K)
    mL = cf_conductor(L)
    mLK = ZZ(mL/mK); # Raises an error if by accident mK \not| mL

    # Construct embeddings K --> L: sigma_s with s mod mLK invertible, s mod mK= 1
    embLK = [ CRT([i,1],[mLK,mK])  for i in range(1,mLK) if gcd(i,mLK)==1 ];

    ext_pol = prod([z - L.gen()**i for i in embLK]);
    ext_pol = Ky( { _key: K(_val) for (_key,_val) in ext_pol.dict().items()} );
    # print("ext_pol", ext_pol);

    # Computation of relative norms by computing the product of conjugates
    xml_pol = ZpX(x**mL-1);
    phiL = ZpX(L.defining_polynomial());
    # [Slow] Compact version, but we need to apply the modulo at each step of the product
    # NA = [ K(list(prod(Zx(_a.polynomial())(x**i).mod(x**mL-1) for i in embLK))[::mLK]) for _a in A];
    
    NA = [ ];
    t = cputime();
    for _a in A:
        _na = ZpX(1);
        _ax = _a.list();
        _ax = [__ax % p**k for __ax in _ax] # do mod p**k
        in_t = cputime();
        for _i in embLK: # Necessary to apply mod x^mL-1 at each step
            _axi = ZpX(sum([_ax[_k]*x**ZZ(mod(_k*_i,mL)) for _k in range(len(_ax))]));
            _na  = (_na*_axi).mod((xml_pol));
        
        in_t = cputime(in_t); # print("inner ", in_t);
        _na = ZpX(_na).mod((phiL));
        _na = Zx(_na)(y).mod(ext_pol);
        assert(_na.degree() == 0);
        
        NA.append(K(_na));
    t = cputime(t);
    print("all relative norms: {:.2f} [{:.2f}/elt]".format(t, t/len(A)));

    # [Test] Resultants --------------------------------------------------------------
    # t = cputime();
    # ZZut.<u,v> = PolynomialRing(ZZ,2);
    # Ext_pol = [ _c.polynomial()(v) for _c in ext_pol.list() ]; assert(len(Ext_pol) == euler_phi(mLK)+1);
    # Ext_pol = ZZut(sum( u**_i*Ext_pol[_i] for _i in range(len(Ext_pol)) ));
    # NA2 = [ ZZut(_a.polynomial()(u)).resultant(Ext_pol,u).mod(K.defining_polynomial()(v)) for _a in A ];
    # NA2 = [ K(_na(1,x)) for _na in NA2 ];
    # t = cputime(t);
    # print("all relative norms using resultants: {:.2f} [{:.2f}/elt]".format(t, t/len(A)));
    # assert(NA == NA2);
    # //--- [END Resultants] Only works for small values, and always twice slower anyway. Does not scale at all.

    return NA;





# ----------------------------------------------------------------------------------
# p-adic lifting (Hensel) -- RNS version


# We compute the "inverse" e-th root: a^(1/e - 1)
# So Let g(x) = a^{e-1} x^e - 1
#        xi+1 = xi - (1/e).xi.g(xi)
# --> Converges to a^{1/e-1}
# --> a.xk = a^{1/e} mod p^k
def cf_padic_lifting_rns(elts, exps, e, p, B, rns_mods):
    r"""
    Compute a p-adic lifting of (prod_i elts[i]^exps[i])^(1/e) in RNS rep.

    INPUT:

    - ``elts`` -- elements of a cyclotomic number field

    - ``exps`` -- integers
    
    - ``e`` -- integer; it is the root exponent

    - ``p`` -- prime integer, inert in elts[i].parent(); modulus for lifting

    - ``B`` -- real number; bound on the norm of the solution

    - ``rns_mods`` -- integers; RNS moduli
    """
    Zx.<x> = PolynomialRing(ZZ);
    K = elts[0].parent();
    m = cf_conductor(K);

    # --------------  WE ASSUME THAT p has been already computed   ------------
    # assert(cf_can_inert(m)); # Needs inert primes (i.e. cyclic ZZ/mZZ*)
    # pp, ee = is_prime_power(m, get_data=True); # Since 1 < m != 2[4], only way is m=4 or pp^ee
    # # An element of order phi(m) ?
    # r = ZZ(IntegerModRing(m).multiplicative_generator());
    # # => inert primes are = r [m] (many choices)
    # _d = 0
    # p = randint(2**60, 2**61)
    # if ee==1: # 1/e does exist mod (p^n - 1)/e
    #     while _d!=1:
    #         p = next_prime_mod_m(m, r, next=p);   print("inert p: ", p, "m: ", m, "r: ", r);
    #         Fp = GF(p, proof=False);
    #         q = ZZ(p**K.degree()); Fq = GF(q, proof=False, modulus=K.polynomial(), name='u'); # if modulus not given, takes forever
    #         _d = gcd(e,(q-1)//e); 
    # else:     # 1/e cannot exist mod (p^n - 1)/e
    #     while _d==0:
    #         p = next_prime_mod_m(m, r, next=p);   print("inert p: ", p, "m: ", m, "r: ", r);
    #         Fp = GF(p, proof=False);
    #         q = ZZ(p**K.degree()); Fq = GF(q, proof=False, modulus=K.polynomial(), name='u'); # if modulus not given, takes forever
    #         _d = gcd(e,(q-1)//e);
    # -------------------------------------------------------------------------
    pp, ee = is_prime_power(mK, get_data=True);
    Fp = GF(p, proof=False);
    q = ZZ(p**K.degree()); Fq = GF(q, proof=False, modulus=K.polynomial(), name='u'); # if modulus not given, takes forever
    print("ff done", flush=True);
    Fpy.<y> = PolynomialRing(Fp); gp = Fpy(K.polynomial());
    Fqz.<z> = PolynomialRing(Fq);
    print("poly done",flush=True);

    # Hensel's height
    log2_k = ceil(max(log(ln(2*B)/ln(p), 2),0)); # [Warn] Protection against small cases giving negative log
    k      = 2**log2_k;
    # print("Hensel's height: {} iterations".format(log2_k));

    Rpk     = IntegerModRing(p**k);
    Rpk.<t> = PolynomialRing(Rpk);
    gt      = Rpk(K.polynomial());

    # Prod mod p^k
    t = cputime();
    se_Pk = prod( power_mod(Rpk(_u.polynomial()), _e, gt) for _u,_e in
                  zip(elts,exps) ).mod(gt);
    a_Pk  = power_mod(se_Pk, e-1, gt); # Say, a = se^{e-1}
    t = cputime(t);
    # print("time se^{{e-1}} mod P^k: {:.2f}".format(t)); 
    # print("se^{{e-1}} mod P^k:", a_Pk);

    # Root in base residue field: first approximation
    t = cputime(); a_P = Fq(Fpy(a_Pk)); t = cputime(t);
    if ee!=1:
        # We should use [AMM77] generalization of Tonnelli-Shanks...
        # Using Cantor-Zassenhaus instead. Seems good enough ?
        pol_pow = a_P*z**e - 1; # [NB] Should compute this mod gt (if e is about not small)
        t = cputime(); is_P  = pol_pow.roots()[0][0]; t = cputime(t); print("roots: {:.2f}".format(t));
    else:
        inv_e = ZZ(-mod(1/e, (q-1)//e)); # careful, this is not mod(-1/e,(q-1)//e) if e not prime
        is_P  = a_P**inv_e;
    # assert(is_P**e * a_P == 1); # Takes time
    # print("a^{{1/e-1}}:", is_P);

    
    # Want to lift root mod P to root mod P^2 etc.

    is_P = Zx(is_P.polynomial());
    a_P = Zx(a_P.polynomial());
    
    _k   = 1; # Just in case log2_k = 0

    for _i in range(log2_k):
        # print("\n"+"-"*15+"Lift stage {}".format(_i)+"-"*15);
        _k = 2**(_i+1);
        Fpky.<yk> = PolynomialRing(IntegerModRing(p**_k));
        gp  = Fpky(gt); # Fpky(K.polynomial()); 
        is_P = Fpky(is_P);
        a_P  = Fpky(a_Pk);

        # print("old sol: {}".format(is_P));
        # xi+1 = xi - 1/e xi f(xi)
        
        
        # # seems to be better to put everything in RNS rep. at each step
        # is_P_rns = [PolynomialRing(IntegerModRing(_p), 'XX')(is_P) for _p in rns_mods]
        # aP_rns = [PolynomialRing(IntegerModRing(_p), 'XX')(a_P) for _p in rns_mods]
        # gt_rns = [PolynomialRing(IntegerModRing(_p), 'XX')(gp) for _p in rns_mods]
        # inv_e_rns = [mod(mod(1/e, p**_k), _p) for _p in rns_mods]
        
        eval_isP = (a_P * power_mod(is_P, e, gp) - 1).mod(gp)
        
        # eval_isP_rns = [(_a_P * power_mod(_is_P, e, _gt) - 1).mod(_gt) for
                        # _a_P,_is_P,_gt in zip(aP_rns, is_P_rns, gt_rns)];
                                
        
        # print("f(xi) mod p^{} = {}".format(_k, eval_isP));

        new_isP = (is_P - mod(1/e, p**_k) * (is_P * eval_isP)).mod(gp);
        # new_isP_rns = [(_is_P - _inv_e * (_is_P * _eval_isP)) for
        #                _is_P,_inv_e,_eval_isP,_gt in
        #                zip(is_P_rns, inv_e_rns, eval_isP_rns, gt_rns)];

        # assert( (power_mod(new_isP, e, gp) * a_P).mod(gp) == 1 ); # Takes time
        is_P = Zx(new_isP);
        # is_P_rns = new_isP_rns;
        
    assert(_k == k);
    
    s_P = (Rpk(is_P)*se_Pk).mod(gt);

    # se_Pk_rns = [PolynomialRing(IntegerModRing(_p), 'XX')(se_Pk) for _p in rns_mods] 
    # s_P_rns = [(_is_P * _se).mod(_gt) for _is_P,_se,_gt in zip(is_P_rns, se_Pk_rns, gt_rns)]

    # print("\n"+"-"*15+"End of p-adic lifting" + "-"*15);
    # assert(power_mod(s_P, e, gt) == se_Pk); # Takes time
    s_P = Zx(s_P);
    s   = K([ (_s if _s < p**k/2. else _s - p**k) for _s in s_P ]);
    # print("SOLUTION: ", s);

    return s;
