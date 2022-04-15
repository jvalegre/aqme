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

calc=crest
flags="$calc"

## order on arguments important for cregen
### Loop over arguments
while [[ $# -gt 0 ]]; do
  case $1 in
      -xyzinbest|--xyzinbest)
        xyzinbest="$2"
        filebest="${xyzinbest%.*}"
        flags="$flags $xyzinbest --cregen"
        shift
        shift
        ;;
      -xyzinall|--xyzinall)
        xyzinall="$2"
        fileall="${xyzinall%.*}"
        flags="$flags $xyzinall"
        shift
        shift
        ;;
    	-ewin|--ewin)
    	  flags="$flags $1 $2"
    	  shift
    	  shift
    	  ;;
      -bthr|--bthr)
    	  flags="$flags $1 $2"
    	  shift
    	  shift
    	  ;;
      -rthr|--rthr)
    	  flags="$flags $1 $2"
    	  shift
    	  shift
    	  ;;
      -ethr|--ethr)
    	  flags="$flags $1 $2"
    	  shift
    	  shift
    	  ;;
    	-t|--time)
    	  taim="$2"
    	  shift
    	  shift
    	  ;;
      -xyzout|--xyzout)
    	  xyzout="$2"
        outfile="${xyzout%.*}".out
    	  shift
    	  shift
    	  ;;
    	-c|--cluster)
    	  flags="$flags $1 $2"
    	  shift
    	  shift
    	  ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
  esac
done

# output
exline="$flags"

echo " "
echo "-excuting $exline"

if [ -z "$fileall" ] && [ -z "$filebest" ]
then
	 echo "NO INPUT!"
else
   $exline > $outfile
   mv crest_ensemble.xyz $xyzout
fi

rm -f scoord* coord* coord.original cregen_0.tmp cregen_1.tmp cre_members *.sorted *.energies struc.xyz wbo gfnff_topo .xcontrol.sample $contoutfile
