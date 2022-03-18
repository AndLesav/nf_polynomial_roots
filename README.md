
# nf_polynomial_roots
Computing polynomial roots in number field through complex embeddings

Code in support of the article "Computing roots of polynomials over number fields using complex embeddings"

https://hal.archives-ouvertes.fr/hal-03608840


## Authors
Andrea Lesavourey.

## Softwares
The code has been tested with:
- Pari/Gp V2.13.1. https://pari.math.u-bordeaux.fr/

- Magma V2.24-10   http://magma.maths.usyd.edu.au/magma/


Plot scripts use Gnuplot 5.4.


## Warning
This code does not come with any guarantee. 


## How to use scripts?

Assume you are in the scripts directory.
You should ensure that directories named "logs", "data", "heads", "inputs" and "figures" are created next to "scripts" before starting.

Bash scripts call on Gp or Magma scripts. Results of computations are printed in the "logs" directory.
If one execute a given bash script  ./Script.sh  without argument, some help will be printed, giving the number of parameters required and a quick description of all parameters (obligatory and optional).


## How to replicate results and plots from the paper?

**CAREFUL** : some parameters may need to be changed ***by hand*** in some plot files

### Heuristics : Section 4
- Figure 1 : Use script `GaussianHeuristic.sh` (computations in Magma) then `plot_gaussian`
- Figures 2 and 3 : Use `SpecLLL.sh` (computations in Magma & Pari/Gp) then `plot_spec_lll`
- Precision from norm, figure 16 : Use `PrecisionEvaluation.sh` then `plot_prec_eval`
- Norm evaluation, figure 4 : Use `NormEvaluation.sh` then `plot_norm_eval`

### Absolute method : Section 5
- Table 1: data are printed for
- Impact of P_K(X), Figure 5: Use `PolfieldField_absolute.sh` then `plot_polfield`
- Impact of s_Z, Figure 6: Use `SizeRoots_absolute.sh` then `plot_sizeroots`
- Impact of #Z_K(f)/deg(f), Figure 7: Use `NbSolutions_absolute.sh` then `plot_nbsol_abs`

### Relative method : Section 5.3
- Impact of heuristics, Figure 8: Use `CompareEqDegree_relative.sh` and `CompareSizeRoots_relative.sh`, then `plot_rel_lll_compar_degeq` and `plot_rel_lll_compar_size`
- Abs vs Rel, #Z(f)/deg(f), Figure 9: Use `CompareNbSol_abs_rel.sh` then `plot_rel_nbsol`
- Abs vs Rel, relative degree [L:K], Figure 10: Use `CompareRelDeg_abs_rel.sh` then `plot_rel_reldeg`
- Cyclotomic fields, Figure 11: Use `Cyclotomics.sh` then `plot_cyclo`
- Kummer fields, Figures 12 and 13: Use `Kummer.sh` then `plot_kummer`


### Summary of the directories from root
 - `./src`: Functions
 - `./scripts`: Bash or Pari/Gp or Magma scripts
 - `./plots`: plot files    
 - `./data`: data obtained
 - `./figures`: where figures plotted from plots are put    
 - `./heads`: directory for temporary files
 - `./logs`: logs
 - `./inputs`: input files computed by Magma scripts
  
License

© 2021, Andrea Lesavourey.
This work is published under the GNU General Public License (GPL) v3.0.
See the LICENSE file for complete statement.
