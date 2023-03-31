#
# This is a template for precision handling
# So far it supports a very simple-minded parametrization:
#
#  .DEFAULT  ---   .4  or .3  (for energies, etc)
#  .FORCES  ---   .4  or .3  (for forces...)
#  .STRESS   ---   .2  or .1  (for stresses, etc)
#
# Atomic section
/ rc = / {printf "%-20s %20.DEFAULTf\n", "rc", $3}
/potential\(screened\)/ {printf "%-20s %20.DEFAULTf\n", "Screened pot", $3}
# basic params
/k-grid: Cutoff / {printf "%-20s %20.DEFAULTf\n", "kgrid-cutoff", $6}
/Internal auxiliary supercell/ {print}
/and projectors:/ {print}
# basic geometry
/Cell vector modules / {printf "%-20s %20.DEFAULTf%20.DEFAULTf%20.DEFAULTf\n", "Cellvec", $7, $8, $9}
/Cell angles /     {printf "%-20s %20.DEFAULTf%20.DEFAULTf%20.DEFAULTf\n", "Cellang", $6, $7, $8}
# energy components
/siesta: Ebs / {printf "%-20s %20.DEFAULTf\n", "Ebs", $4}
/siesta: Eions /  {printf "%-20s %20.DEFAULTf\n", "Eions", $4}
/siesta: Ena   /  {printf "%-20s %20.DEFAULTf\n", "Ena", $4}
/siesta: Ekin  /  {printf "%-20s %20.DEFAULTf\n", "Ekin", $4}
/siesta: Enl   / {printf "%-20s %20.DEFAULTf\n", "Enl", $4}
/siesta: Eharris / {printf "%-20s %20.DEFAULTf\n", "Eharris", $4}
/siesta: Etot / {printf "%-20s %20.DEFAULTf\n", "Etot", $4}
/siesta: FreeEng / {printf "%-20s %20.DEFAULTf\n", "FreeEng", $4}
/siesta: E_KS\(eV\)/  {printf "%-20s %20.DEFAULTf\n", "KS E", $4}
/siesta: Eharris\(eV\)/  {printf "%-20s %20.DEFAULTf\n", "Harris E", $4}
#forces
/sqrt\( Sum f/ {printf "%-20s %20.FORCESf\n", "Res force", $2}
/constrained/ {printf "%-20s %20.FORCESf\n", "Max force", $2}
#pressure
/Stress-tensor-Voigt/ {printf "%-20s %10.STRESSf%10.STRESSf%10.STRESSf%10.STRESSf%10.STRESSf%10.STRESSf\n", "Stress", $3, $4, $5, $6, $7, $8}
# polarization
/Along the lattice vectors / {printf "%-20s %20.DEFAULTf%20.DEFAULTf%20.DEFAULTf\n", "Pol-latt", $6, $7, $8}
/Along cartesian / {printf "%-20s %20.DEFAULTf%20.DEFAULTf%20.DEFAULTf\n", "Pol-cart", $5, $6, $7}
# optical
/Checking f-sum rule/  {printf "%-20s %20.DEFAULTf\n", "f-sum rule", $5}
# spin
/\(Qup-Qdown\)/ {printf "%-20s %20.DEFAULTf\n", "Spin-pol", $7}
# mulliken
/mulliken: Qtot / {printf "%-20s %20.DEFAULTf\n", "Spin-pol", $4}




