#!/bin/bash

# default variables - machine dependent, uncomment as appropriate!
runcrest=/usr/local/xtb/crest

charge=${charge:-0}
xyzoutall=${xyzoutall:-crest_conformers.xyz}
xyzoutbest=${xyzoutbest:-crest_best.xyz}



# input geometry (xyz)
input=$1
file="${input%.*}"

echo -e "-  RUNNING $file WITH crest \c"

# output
outfile="${xyzoutall%.*}".out
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

# run xTB if xyz file supplied
if [ -z "$file" ]
then
	 echo "NO INPUT!"
else
	 $runcrest $file.xyz -c $charge -ewin $ewin > $outfile && mv crest_conformers.xyz $xyzoutall && mv crest_best.xyz $xyzoutbest
fi

rm -f scoord* coord* coord.original cregen_0.tmp cregen_1.tmp *.sorted cre_members crest.energies crest_* struc.xyz wbo gfnff_topo .xcontrol.sample $contoutfile scoord*
