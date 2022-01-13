%nprocshared=8
%mem=16GB
# pop=(nbo6read,savenbos) wb97xd/gen scrf=(smd,solvent=dichloromethane)

H2

0 1
 H   0.00000000   0.00000000   0.37011200
 H   0.00000000   0.00000000  -0.37011200

H 0
6-31G*
****

$nbo bndidx $end

