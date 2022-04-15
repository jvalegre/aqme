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

xyzoutall=${xyzoutall:-crest_conformers.xyz}
xyzoutbest=${xyzoutbest:-crest_best.xyz}

calc=crest
flags="$calc"

### Loop over arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input)
      file="$2"
	  flags="$flags $file"
      shift # past argument
      shift # past value
      ;;
    	-c|--charge)
    	  charge="$2"
    	  flags="$flags -chrg $charge"
    	  shift
    	  shift
    	  ;;
    	-u|--uhf)
    	  mult="$2"
    	  flags="$flags -uhf $mult"
    	  shift
    	  shift
    	  ;;
    	-s|--solvent)
    	  solvent="$2"
    	  flags="$flags --alpb $solvent"
    	  shift
    	  shift
    	  ;;
    	-p|--nproc)
    	  proc="$2"
    	  flags="$flags -T $proc"
    	  shift
    	  shift
    	  ;;
    	-fz|--constrain)
    	  constrain="$2"
    	  rm -f constrain.inp .xcontrol.sample
    	  crest coord --constrain $constrain > const.out
    	  sed -e 's/ > //g' const.out | tail -8 > constrain.inp
    	  flags="$flags -I constrain.inp"
    	  rm -f const.out
    	  shift
    	  shift
    	  ;;
    	-ewin|--ewin)
    	  flags="$flags $1 $2"
    	  shift
    	  shift
    	  ;;
      -cbonds|--cbonds)
    	  flags="$flags $1 $2"
    	  shift
    	  shift
    	  ;;
    	--nci)
    	  flags="$flags $1"
    	  shift
    	  ;;
    	--noreftopo)
    	  flags="$flags $1"
    	  shift
    	  ;;
    	--protonate)
    	  flags="$flags $1"
    	  shift
    	  ;;
    	--deprotonate)
    	  flags="$flags $1"
    	  shift
    	  ;;
    	--tautomerize)
    	  flags="$flags $1"
    	  shift
    	  ;;
    	-t|--time)
    	  taim="$2"
    	  shift
    	  shift
    	  ;;
      -xyzoutall|--xyzoutall)
    	  xyzoutall="$2"
        outfile="${xyzoutall%.*}".out
    	  shift
    	  shift
    	  ;;
      -xyzoutbest|--xyzoutbest)
    	  xyzoutbest="$2"
    	  shift
    	  shift
    	  ;;
    	-a|--add)
    	  add="$2"
    	  flags="$flags $add"
    	  shift
    	  shift
    	  ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
  esac
done

### If constraints were set, convert coordinate file format to Turbomole
if [ ! -z "$constrain"] && [ "$file" =! "coord" ]; then
    t2x $file > coord
	flags=$(echo $flags | sed "s/$file/coord/")
fi

# output
exline="$flags"

echo " "
echo "-excuting $exline"

if [ -z "$file" ]
then
	 echo "NO INPUT!"
else
   $exline > $outfile
   mv crest_conformers.xyz $xyzoutall && mv crest_best.xyz $xyzoutbest
fi

rm -f coord* coord.original cregen_0.tmp cregen_1.tmp *.sorted cre_members crest.energies struc.xyz wbo gfnff_topo .xcontrol.sample $contoutfile scoord*
