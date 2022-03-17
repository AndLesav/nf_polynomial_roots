# nf_polynomial_roots
Computing polynomial roots in number field through complex embeddings

Code in support of the article "Computing roots of polynomials over number fields using complex embeddings

"https://hal.archives-ouvertes.fr/hal-03608840


## Authors
Andrea Lesavourey.

## Softwares
The code has been tested with:
    - Pari/Gp V2.13.1. https://pari.math.u-bordeaux.fr/
    - Magma V2.24-10   http://magma.maths.usyd.edu.au/magma/


## Warning
This code does not come with any guarantee. 


## How to use scripts?

Assume you are in the scripts directory.
You should ensure that directories named "logs", "data", "heads", "inputs" and "figures" are created next to "scripts" before starting.

Bash scripts call on Gp or Magma scripts. Results of computations are printed in the "logs" directory.
If one execute a given bash script  ./Script.sh  without argument, some help will be printed, giving the number of parameters required and a quick description of all parameters (obligatory and optional).




### Summary of the directories from root
    "./src": Functions
    "./scripts": Bash or Pari/Gp or Magma scripts
    "./plots": plot files    
    "./data": data obtained
    "./figures": where figures plotted from plots are put    
    "./heads": directory for temporary files
    "./logs": logs
    "./inputs": input files computed by Magma scripts
  



License

© 2021, Andrea Lesavourey.
This work is published under the GNU General Public License (GPL) v3.0.
See the LICENSE file for complete statement.
