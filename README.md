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


***WARNING : description does not match new versions of the paper or the code***

***TODO : modify the description accordingly***


Assume you are in the scripts directory.
You should ensure that directories named "logs", "data", "heads", "inputs" and 
"figures" are created next to "scripts" before starting.

Bash scripts call on Gp or Magma scripts. Results of computations are printed 
in the "logs" directory.

If one execute a given bash script  ./Script.sh  without argument, some help 
will be printed, giving the number of parameters required and a quick 
description of all parameters (obligatory and optional).


## How to replicate results and plots from the paper?

**CAREFUL** : some parameters may need to be changed ***by hand*** in some plot files

### Heuristics:

#### Section 3:
- Gaussian heuristic, figure 1: use script `GaussianHeuristic.sh` (computations in Magma) then
  `plot_gaussian`

#### Section 5:
- SpecLLL, figures 2 and 3: use `SpecLLL.sh` (computations in Magma & Pari/Gp) then
`plot_spec_lll`

- Precision from norm, figure 16 (appendix): use `PrecisionEvaluation.sh` then
`plot_prec_eval`

- Norm evaluation, figure 4: use `NormEvaluation.sh` then `plot_norm_eval`

- Early abort in decoding, table 1: use `EarlyAbort.sh` 


### "Good" and "bad" fields : Section 6.1.
- Statistical data, table 2: from previous version of `Polfield_absolute.sh`

- Proportion of bad fields, table 3: use `PropBadFields.sh`

- Timings for good / bad fields, table 4: use `TimingsBadFields.sh`


### Absolute method : Section 6.2.
- Impact of $P_K(X)$, figure 5: Use `PolfieldField_absolute.sh` then 
`plot_polfield`

- Impact of $s_Z$, figure 6: Use `SizeRoots_absolute.sh` then `plot_sizeroots`



### Relative method : Section 5.3

Some scripts comparing absolute and relative methods require to pre-compute 
defining polynomials which are then saved in the `./inputs` directory.
This can be done using the scripts whose names are starting by `Field_creation`.

- Impact of heuristics, figure 7: use `CompareEqDegree_relative.sh` and
`CompareSizeRoots_relative.sh`, then `plot_rel_lll_compar_degeq` and 
`plot_rel_lll_compar_size`

- Abs vs Rel, $|Z(f)|/\deg(f)$, figure 8: use `FieldCreation_nbsol.sh` to 
	create number fields, then `CompareNbSol_abs_rel.sh` and `plot_rel_nbsol`

- Abs vs Rel, relative degree $[L:K]$, figure 9: use `FieldCreation_reldeg.sh`
to create number fields, then `CompareRelDeg_abs_rel.sh` and `plot_rel_reldeg`

- Cyclotomic fields with extension $K_m/K_m^+$, figure 10: use `Cyclotomics.sh`
then `plot_cyclo`

- Prime powers cyclotomic fields, table 5 : use `PrimePowerCyclotomics.sh`

- Kummer fields, figures 11 and 12: use `FieldCreation_kummer.sh` to create 
fields, then `Kummer.sh` (with version split and single respectively) and
`plot_kummer_general`

- Small degree number fields and large roots, figure 13: use 
`FieldCreation_reldeg_fixed.sh` to create fields, then `RelativeDegreeFixed.sh` 
and `plot_reldeg_fixed`


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
