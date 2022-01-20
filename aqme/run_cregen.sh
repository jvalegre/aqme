#!/bin/bash

# default variables - machine dependent, uncomment as appropriate!
runcrest=/usr/local/xtb/crest

charge=${charge:-0}
xyzout=${xyzout:-crest_conformers.xyz}
ethr=${ethr:-0.2}
rthr=${rthr:-0.125}
bthr=${bthr:-0.01}
ewin=${ewin:-6}


# input geometry (xyz)
inputall=$1
inputbest=$2
fileall="${inputall%.*}"
filebest="${inputbest%.*}"
echo -e "-  RUNNING $fileall and $filebest WITH cregen \c"

# output
outfile="${xyzout%.*}".out
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
if [ -z "$fileall" ] && [ -z "$filebest" ]
then
	 echo "NO INPUT!"
else
   $runcrest $filebest.xyz -cregen $fileall.xyz -ethr $ethr -rthr $rthr -bthr $bthr -ewin $ewin > $outfile  && mv crest_ensemble.xyz $xyzout
fi

rm -f scoord* coord* coord.original cregen_0.tmp cregen_1.tmp cre_members *.sorted *.energies crest_* struc.xyz wbo gfnff_topo .xcontrol.sample $contoutfile
