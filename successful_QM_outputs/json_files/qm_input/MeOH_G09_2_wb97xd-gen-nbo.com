%nprocshared=8
%mem=16GB
# pop=(nbo6read,savenbos) wb97xd/gen scrf=(smd,solvent=dichloromethane)

MeOH_G09_2

0 1
 C   0.65978700  -0.01959100   0.00000000
 H   1.08588800   0.98380100   0.00000100
 H   1.02538400  -0.54453500  -0.88940900
 H   1.02538400  -0.54453700   0.88940800
 O  -0.74538700   0.12221100   0.00000000
 H  -1.13228700  -0.75487100   0.00000000

H O 0
6-31G*
****
C 0
def2svp
****

$nbo bndidx $end

