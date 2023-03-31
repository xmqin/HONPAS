#
# Simple awk script to create a "digest" of key results in a
# SIESTA .out file.
# This digest can then be compared with a reference to help in
# finding discrepancies.
# The tolerances (widths of decimal fields) in this version are
# hard-wired.
#
# A. Garcia, Feb. 2011
#
# Atomic section
#
/ rc = / {printf "%-20s %20.4f\n", "rc", $3}
/potential\(screened\)/ {printf "%-20s %20.4f\n", "Screened pot", $3}
#
# basic params
#
/k-grid: Cutoff / {printf "%-20s %20.4f\n", "kgrid-cutoff", $6}
/Internal auxiliary supercell/ {print}
/and projectors:/ {print}
#
# basic geometry
#
/Cell vector modules / {printf "%-20s %20.4f%20.4f%20.4f\n", "Cellvec", $7, $8, $9}
/Cell angles /     {printf "%-20s %20.4f%20.4f%20.4f\n", "Cellang", $6, $7, $8}
#
# energy components
#
/siesta: Ebs / {printf "%-20s %20.4f\n", "Ebs", $4}
/siesta: Eions /  {printf "%-20s %20.4f\n", "Eions", $4}
/siesta: Ena   /  {printf "%-20s %20.4f\n", "Ena", $4}
/siesta: Ekin  /  {printf "%-20s %20.4f\n", "Ekin", $4}
/siesta: Edftu /  {printf "%-20s %20.4f\n", "Edftu", $4}
/siesta: Eso   /  {printf "%-20s %20.4f\n", "Eso", $4}
/siesta: Enl   / {printf "%-20s %20.4f\n", "Enl", $4}
/siesta: Eharris / {printf "%-20s %20.4f\n", "Eharris", $4}
/siesta: Etot / {printf "%-20s %20.4f\n", "Etot", $4}
/siesta: FreeEng / {printf "%-20s %20.4f\n", "FreeEng", $4}
/siesta: E_KS\(eV\)/  {printf "%-20s %20.4f\n", "KS E", $4}
/siesta: Eharris\(eV\)/  {printf "%-20s %20.4f\n", "Harris E", $4}
#
#forces
#
/sqrt\( Sum f/ {printf "%-20s %20.4f\n", "Res force", $2}
/constrained/ {printf "%-20s %20.4f\n", "Max force", $2}
#
#pressure
#
/Stress-tensor-Voigt/ {printf "%-20s %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n", "Stress", $3, $4, $5, $6, $7, $8}
#
# polarization
#
/Along the lattice vectors / {printf "%-20s %20.4f%20.4f%20.4f\n", "Pol-latt", $6, $7, $8}
/Along cartesian / {printf "%-20s %20.4f%20.4f%20.4f\n", "Pol-cart", $5, $6, $7}
#
# optical
#
/Checking f-sum rule/  {printf "%-20s %20.4f\n", "f-sum rule", $5}
#
# spin
#
/\(Qup-Qdown\)/ {printf "%-20s %20.4f\n", "Spin-pol", $7}
#
# mulliken
#
/mulliken: Qtot / {printf "%-20s %20.4f\n", "Spin-pol", $4}
