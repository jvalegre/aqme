#!/bin/bash

# default variables - machine dependent, uncomment as appropriate!
# runcrest=/usr/local/xtb/crest
XTBHOME=/usr/local/xtb
export XTBPATH=${XTBHOME}/share/xtb:${XTBHOME}:${HOME}
export PATH=${PATH}:${XTBHOME}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${XTBHOME}/lib
export PYTHONPATH=${PYTHONPATH}:${XTBHOME}/python
export PATH=/usr/local/xtb:$PATH

ulimit -s unlimited
export OMP_STACKSIZE=1G
export OMP_NUM_THREADS=12,1
export OMP_MAX_ACTIVE_LEVELS=1
export MKL_NUM_THREADS=12

charge=${charge:-0}
nproc=${nproc:-24}
force_const=${force_const:-1.0}
max_cycle=${max_cycle:-100}
json=${json:-true}
xyzoutall=${xyzoutall:-crest_conformers.xyz}
xyzoutbest=${xyzoutbest:-crest_best.xyz}
constraint=${constraint:-}
ewin=${ewin:-6}
ethr=${ethr:-0.2}
rthr=${rthr:-0.125}
bthr=${bthr:-0.01}
bthr=${bthr:-0.02}

# input geometry (xyz)
input=$1
file="${input%.*}"

echo -e "-  RUNNING $file WITH crest \c"

# output
outfile="${xyzoutall%.*}".out
cregenoutfile="${xyzoutall%.*}".cregen.out
contoutfile="${xyzoutall%.*}".constraint.out
echo -e "-  outfile $outfile \c"

# user-defined variables
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        echo -e $1 $2 "\c"
   fi
  shift
done

echo -e "-  constraint $constraint \c"
# run xTB if xyz file supplied
if [ -z "$file" ]
then
	 echo "NO INPUT!"
else
	 crest $file.xyz --constrain $constraint > $contoutfile
	 crest $file.xyz --cinp .xcontrol.sample -c $charge -ewin $ewin --noreftopo --nci --cbonds $cbonds > $outfile && mv crest_conformers.xyz $xyzoutall && mv crest_best.xyz $xyzoutbest
fi

rm -f coord* coord.original cregen_0.tmp cregen_1.tmp *.sorted cre_members crest.energies crest_* struc.xyz wbo gfnff_topo .xcontrol.sample $contoutfile scoord*
