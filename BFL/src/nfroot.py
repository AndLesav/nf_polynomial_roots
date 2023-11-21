# Code in support of ePrint:2021/1384
# Copyright 2021, Andrea Lesavourey
# GPL-3.0-only (see LICENSE file)


import sys
TWSTI_PATH="../Tw-Sti/"; 
sys.path.append(TWSTI_PATH + "src/");

from sage.all import *
import fp
from pideal import *
from lattice import *
from nf import *
proof.number_field(False);

# sort family of nf elements according to coefficient norm
# Hv = elements in coeff representation
def __nf_coeff_sort(Hv, K):
    sort_idx = sorted(range(len(Hv)), key=lambda idx: Hv[idx].norm());
    Hv.sort(key=norm);
    return [K(list(h)) for h in Hv], sort_idx;

# sort family of nf elements according to coefficient norm
# return a list of index and the log-norms of the basis
def __nf_coeff_sort_factored(B, exps):
    Bn = [log(vector(_b).norm(), 2) for _b in B] # log norms of basis
    print (Bn)
    # norm of elts corresponding to (B,exps)
    Sn = [sum(abs(_e) * _n for _e,_n in zip(_exp, Bn)) for _exp in exps]
    sort_idx = sorted(range(len(Sn)), key=lambda idx: Sn[idx]); # sort
    return sort_idx, Sn

# evaluation of precision needed
def __prec_eval_sqrtn(g, n):
    return len(g)*ceil(log(g.norm(), 2))//n;


# evaluation of precision needed -- factored form
def __prec_eval_sqrtn_factored(dim, bn, n):
    return dim * ceil(bn) // n;

# evaluation of precision needed -- on a family
def __prec_eval_sqrtn_fam(G, n):
    return max([len(g)*ceil(log(g.norm(), 2)) // n for g in G]);


def row_augment(M, N):
    return (M.transpose().augment(N.transpose())).transpose()

# special lll for normal bases (power basis ?)
def spec_lll(M):
    dim = M.nrows(); coeff = dim // 2;
    r = M.ncols() - M.nrows();
    
    N = copy(M[:2,:r+2]);
    N[1, 1 + r] *= N[0].norm()**2;
    N = N.LLL();
    
    for i in range(2, dim):     # we size-reduce M[i-1] wrt M[:i-2], now we need to add M[i]
        N[i-1, i-1+r] = coeff;
        N = row_augment(N, M[i,:N.ncols()]);
        N = N.augment(zero_matrix(ZZ, i+1, 1));  N[i, i+r] = coeff;
        N[i,:r] *= coeff;
        for j in range(i-1):
            N[i,:r] += N[i-1,r+j] * M[j+1,:r]
            N[i, r+j+1] = N[i-1, r+j]
        N[i,:r] /= coeff
        N[i, i+r] *= N[i].norm()**2
        N = N.LLL()
    
    N[dim-1, r+ dim-1] = coeff;
    return N.LLL();
        

# creation of reduced basis for given precision
def reduced_basis_latt(K, Binf, bprec):
    dim = K.degree(); 
    B = [-round( (2**bprec*Binf[i]).real_part()) for i in range(dim//2)];
    coeff = dim // 2;
    M = coeff*identity_matrix(dim//2, ZZ);
    M = M.insert_row(0, B)
    M = M.transpose()
    # return lll_ZZ(M)[0];
    # return M.LLL()
    return spec_lll(M);

def init_basis_latt(K):
    dim = K.degree(); 
    B = vector(ZZ, [-1 for i in range(dim//2)]);
    coeff = dim // 2;
    
    M = coeff*identity_matrix(dim//2, ZZ);
    M = M.insert_row(0, B)
    M = M.transpose()
    return M;


# use previous unimodular matrix U to start the reduction
def update_basis_latt(K, U, Binf, bprec):
    dim = K.degree(); 
    B = [-round(((2**bprec)*(Binf[i].real_part()))) for i in range(dim//2)];
    coeff = dim // 2;
    V = U/coeff;
    M = coeff*identity_matrix(dim//2, ZZ);
    M = M.insert_row(0, B)
    M = M.transpose()
    # return spec_lll(V*M);
    # return lll_ZZ(V*M)[0];
    if (V == identity_matrix(ZZ, dim)): return spec_lll(M);
    return (V*M).LLL();

# decoding by nearest plane and kannan's embedding technique
# for better performances and accuracy replace by proper NP
def test_decode(BL, x_emb, bprec):
    dim = BL.nrows()
    coeff = dim 
    r = BL.ncols()-BL.nrows()

    Y = [round((2**(bprec))*x_emb[i]) for i in range(r)]; # embeddings of target
    Z = list(zero_vector(ZZ, dim))
    Y += Z
    
    M = (BL.transpose()).augment(vector(ZZ, Y)).transpose()
    m = max(map(abs, M.coefficients()))
    Z += [2**(dim//2) * m ];    # large coeff for Kannan's technique
    
    M = M.augment(vector(ZZ,Z));

    # M , _ = lll_ZZ(M);
    M = M.LLL()
    return M[dim][r:M.ncols()-1]/coeff


# d-th root computation of h in K
# assume L=basis lattice for decoding is precomputed
# assume M=inv relative mink precomputed as well (2x2 matrix)
def cf_sqrtn(h, d, Minf, K, L, pinf, Binf, bprec):
    zeta = K.gen();
    eta = zeta + 1/zeta;
    hinf = pinf[0](h); # put in CC
    hconj = hinf.conjugate()
    ginf = hinf**(1/d);
    gconj = ginf.conjugate();

    bprec_max = pinf[0].codomain().precision()
    zeta_d = CyclotomicField(d).gen().complex_embedding(bprec_max)    
    bprec_new = bprec;
    b = False;

    while (not b):
        
        if bprec_new > bprec_max:
            # print("updating pinf and all subsequent objects");
            bprec_max = 2*bprec_new;
            pinf = get_inf_places(K, bprec_max);
            Binf = [pinf[0](eta**k) for k in range(K.degree()//2)];
            zinf = pinf[0](zeta);
            zeta_d = CyclotomicField(d).gen().complex_embedding(bprec_max)
            Minf = Matrix([[1, zinf], [1, zinf**(-1)]])**(-1);
            hinf = pinf[0](h); # put in CC
            ginf = hinf**(1/d);

        if bprec_new != bprec:
            print("\tin nfroots: prec update [{} --> {}]".format(bprec, bprec_new));
            U = L[:,1:];
            s = cputime()
            L = update_basis_latt(K, U, Binf, bprec_new).change_ring(ZZ);
            bprec = bprec_new;

        i = 0
        while (not b) and (i < d):
           
            i += 1
            ginf *= zeta_d
            gconj = ginf.conjugate();

            # going to rep of g as g = g0 + g1 zeta_m
            ginf_1  = vector([ginf, gconj]);
            ginf_1 = Minf*(ginf_1.column()); 

            # determine coeffs of g0 wrt the basis of K^+ 
            g0 = test_decode(L, [ginf_1[0,0].real_part()], bprec);
            g0 = sum([g0[k]*eta**k for k in range(len(g0))]);

            # determine coeffs of g1 wrt the basis of K^+ 
            g1 = test_decode(L, [ginf_1[1,0].real_part()], bprec);
            g1 = sum([g1[k]*eta**k for k in range(len(g1))]);

            # modify this ? maybe too long
            b = ((g0+g1*zeta)**d==h); 
            print("\tin nfroots: one root tested (try {})".format(i));
            
        bprec_new = ceil(sqrt(2)*bprec); # update precision if needed
        
    return g0+g1*zeta, L, bprec, pinf, Binf, Minf;

# d-th root computation of family in K
# [Nb] Now the order of the output is guaranteed to be the same as the input
def cf_sqrtn_fam(H, d, K):
    dim = K.degree();
    zeta = K.gen();
    eta = zeta + 1/zeta;
    pinf = get_inf_places(K)    # initial values for embeddings
    
    # first embeddings to determine right precision
    zinf = pinf[0](zeta);
    Minf = Matrix([[1, zinf], [1, zinf**(-1)]])**(-1);
    N = [vector(eta**k) for k in range(dim//2)] + [vector(zeta*eta**k) for k in range(dim//2)];
    N = Matrix(N);

    # sort elements by coeff norm
    Hv = [vector(H[i]) for i in range(len(H))];
    Hn, sort_idx = __nf_coeff_sort(Hv, K);
    # evaluation of max precision needed : ad-hoc for now
    bprec_max = __prec_eval_sqrtn_fam(Hv, d) + 2 * ceil(log(Minf.norm(), 2) * dim*log(dim) + log(dim) * log(log(dim)) * dim + log(dim) * dim + log(N.norm(), 2) );
    pinf = get_inf_places(K, bprec_max)
    
    # embeddings of bases that we will use
    zinf = pinf[0](zeta);
    Minf = Matrix([[1, zinf], [1, zinf**(-1)]])**(-1);
    Binf = [pinf[0](eta**k) for k in range(K.degree()//2)];

    L = init_basis_latt(K);
    R = [];

    bprec = 0;
    
    for k in range(len(Hn)):

        # precision need for current element : a bit ad-hoc as well 
        bprec_new = __prec_eval_sqrtn(vector(Hn[k]), d) // 2 + ceil(log(N.norm(), 2)*dim // 3 + log(dim) * dim);
        
        # update basis lattice that we use if necessary
        if bprec_new > bprec:
            print("Update prec to {}".format(bprec_new), end='', flush=True);
            bprec = bprec_new;
            t = walltime();
            L = update_basis_latt(K, L[:,1:], Binf, bprec);
            t = walltime(t);
            print("\t[done] t={:.2f}".format(t), flush=True);
            
        print("Computing root(#{})\n\tlog_2 red norm {:.2f}\n".format(k, float(log(t2_norm(Hn[k]),2))), end='', flush=True);
        t = walltime()
        r, L, bprec, pinf, Binf, Minf = cf_sqrtn(Hn[k], d, Minf, K, L, pinf, Binf, bprec);
        t = walltime(t);
        print("\x1b[32m[OK]\x1b[0m {:.2f} t={:.2f}\n".format(float(log(t2_norm(r),2)),t), end='', flush=True);
        assert(r**d == Hn[k]);
        R += [r]

    _, R = zip(*sorted(zip(sort_idx, R))); # Reordering
    R    = list(R);
    # assert(all(R[k]**d == H[k] for k in range(len(H)))); # Trust but control

    return(R);

# #############################################
# try for factored form

# d-th root computation of h in K
# assume L=basis lattice for decoding is precomputed
# assume M=inv relative mink precomputed as well (2x2 matrix)
def cf_sqrtn_factored(B, exp, d, Minf, K, L, pinf, Binf, bprec):
    r"""
    Computes the d-th root of prod_i B[i]^exp[i] in K
    
    INPUT :
    
    - ``B`` -- factor basis
    - ``exp`` -- list of exponents corresponding to the decomposition of s in the factor basis
    - ``d`` -- exponent for the root computation, i.e., s = r**d for some r

    /!\ Careful : here two standard of notation collide
    ``Binf`` -- basis of K^+ put into CC by the action of pinf[0] ; it is NOT the embeddings
    of the elements of B (which is the factor basis)
    """
    zeta = K.gen();
    eta = zeta + 1/zeta;
    C = [pinf[0](_b) for _b in B]
    hinf = prod(_c**_e for _c,_e in zip(C,exp)) # put s in CC
    hconj = hinf.conjugate()
    ginf = hinf**(1/d);
    gconj = ginf.conjugate();

    bprec_max = pinf[0].codomain().precision()
    zeta_d = CyclotomicField(d).gen().complex_embedding(bprec_max)    
    bprec_new = bprec;
    b = False;

    while (not b):
        
        if bprec_new > bprec_max:
            # print("updating pinf and all subsequent objects");
            bprec_max = 2*bprec_new;
            pinf = get_inf_places(K, bprec_max);
            Binf = [pinf[0](eta**k) for k in range(K.degree()//2)];
            zinf = pinf[0](zeta);
            zeta_d = CyclotomicField(d).gen().complex_embedding(bprec_max)
            Minf = Matrix([[1, zinf], [1, zinf**(-1)]])**(-1);
            C = [pinf[0](_b) for _b in B]
            hinf = prod(_c**_e for _c,_e in zip(C,exp)) # put s in CC
            ginf = hinf**(1/d);

        if bprec_new != bprec:
            print("\tin nfroots: prec update [{} --> {}]".format(bprec, bprec_new));
            U = L[:,1:];
            s = cputime()
            L = update_basis_latt(K, U, Binf, bprec_new).change_ring(ZZ);
            bprec = bprec_new;

        i = 0
        while (not b) and (i < d):
           
            i += 1
            ginf *= zeta_d
            gconj = ginf.conjugate();

            # going to rep of g as g = g0 + g1 zeta_m
            ginf_1  = vector([ginf, gconj]);
            ginf_1 = Minf*(ginf_1.column()); 

            # determine coeffs of g0 wrt the basis of K^+ 
            g0 = test_decode(L, [ginf_1[0,0].real_part()], bprec);
            g0 = sum([g0[k]*eta**k for k in range(len(g0))]);

            # determine coeffs of g1 wrt the basis of K^+ 
            g1 = test_decode(L, [ginf_1[1,0].real_part()], bprec);
            g1 = sum([g1[k]*eta**k for k in range(len(g1))]);

            # ----- check if we got a solution ---------------------------------
            # modify this ? 
            b = (pinf[0](g0+g1*zeta) -  hinf) < 10**(-bprec // 5)
            # ------------------------------------------------------------------

            print("\tin nfroots: one root tested (try {})".format(i));
            
        bprec_new = ceil(sqrt(2)*bprec); # update precision if needed
        
    return g0+g1*zeta, L, bprec, pinf, Binf, Minf;


# d-th root computation of family in K
# [Nb] Now the order of the output is guaranteed to be the same as the input
def cf_sqrtn_fam_factored(B, exps, d, K):
    dim = K.degree();
    zeta = K.gen();
    eta = zeta + 1/zeta;
    pinf = get_inf_places(K)    # initial values for embeddings
    
    # -----   first embeddings to determine right precision ------
    zinf = pinf[0](zeta);
    Minf = Matrix([[1, zinf], [1, zinf**(-1)]])**(-1);
    N = [vector(eta**k) for k in range(dim//2)] + [vector(zeta*eta**k) for k in range(dim//2)];
    N = Matrix(N);
    # ------------------------------------------------------------
    
    # -----  sort elements by coeff norm and compute max precision to use -----
    sort_idx, Bn = __nf_coeff_sort_factored(B, exps) # Bn = log-norms of basis elts

    print (Bn)
    
    # evaluation of max precision needed : ad-hoc for now
    bprec_max =  2 * __prec_eval_sqrtn_factored(dim, Bn[sort_idx[-1]], d) + 2 * ceil(log(Minf.norm(), 2) * dim*log(dim) + log(dim) * log(log(dim)) * dim + log(dim) * dim + log(N.norm(), 2) );
    pinf = get_inf_places(K, bprec_max)
    # -------------------------------------------------------------------------

    # ----- embeddings of bases that we will use and initialisation ------------
    zinf = pinf[0](zeta);
    Minf = Matrix([[1, zinf], [1, zinf**(-1)]])**(-1);
    Binf = [pinf[0](eta**k) for k in range(K.degree()//2)];
    L = init_basis_latt(K);
    R = [];
    bprec = 0;
    # --------------------------------------------------------------------------

    
    for k in range(len(exps)):

        # precision need for current element : a bit ad-hoc as well 
        bprec_new = 2 * __prec_eval_sqrtn_factored(dim, Bn[sort_idx[k]], d) // 2 + ceil(log(N.norm(), 2)*dim // 3 + log(dim) * dim);
        
        # update basis lattice that we use if necessary
        if bprec_new > bprec:
            print("Update prec to {}".format(bprec_new), end='', flush=True);
            bprec = bprec_new;
            t = walltime();
            L = update_basis_latt(K, L[:,1:], Binf, bprec);
            t = walltime(t);
            print("\t[done] t={:.2f}".format(t), flush=True);
            
        # print("Computing root(#{})\n\tlog_2 red norm {:.2f}\n".format(k, float(log(t2_norm(Hn[k]),2))), end='', flush=True);
        t = walltime()
        r, L, bprec, pinf, Binf, Minf = cf_sqrtn_factored(B, exps[sort_idx[k]], d, Minf, K, L, pinf, Binf, bprec);
        t = walltime(t);
        print("\x1b[32m[OK]\x1b[0m {:.2f} t={:.2f}\n".format(float(log(t2_norm(r),2)),t), end='', flush=True);
        # assert(r**d == Hn[k]);
        R += [r]

    _, R = zip(*sorted(zip(sort_idx, R))); # Reordering
    R    = list(R);
    # assert(all(R[k]**d == H[k] for k in range(len(H)))); # Trust but control

    return(R);
