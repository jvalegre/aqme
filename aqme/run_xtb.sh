#!/bin/bash

# default variables - machine dependent, uncomment as appropriate!
# runxtb=/usr/local/xtb/bin/xtb
#runxtb=/usr/local/xtb_exe/bin/xtb
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

xyzout=${xyzout:-xtbopt.xyz}
force_const=${force_const:-1.5}
max_cycle=${max_cycle:-300}

echo -e "\$constrain" >> constrain.inp

calc=xtb
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
	-j|--job)
      job="$2"
	  flags="$flags --$job"
      shift # past argument
      shift # past value
      ;;
	-c|--charge)
	  charge="$2"
	  flags="$flags -c $charge"
	  shift
	  shift
	  ;;
	-u|--uhf)
	  mult="$2"
	  flags="$flags -u $mult"
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
	  flags="$flags -P $proc"
	  shift
	  shift
	  ;;
  --dist)
    dist="$2"
    echo -e "distance: $dist" >> constrain.inp
    shift
	  shift
	  ;;
  --angle)
    angle="$2"
    echo -e "angle: $angle" >> constrain.inp
    shift
	  shift
	  ;;
  --dihedral)
    dihedral="$2"
    echo -e "dihedral: $dihedral" >> constrain.inp
    shift
	  shift
	  ;;

  -xyzout|--xyzout)
    xyzout="$2"
    outfile="${xyzout%.*}".out
    shift
    shift
    ;;

	-t|--time)
	  taim="$2"
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

# write constraints to file
echo -e "force constant=$force_const" >> constrain.inp
echo -e "\$opt\nmaxcycle=$max_cycle\n\$end" >> constrain.inp

flags="$flags -I constrain.inp"

# output
exline="$flags"

echo " "
echo "-excuting $exline"

if [ -z "$file" ]
then
	 echo "NO INPUT!"
else
   $exline > $outfile && mv xtbopt.xyz $xyzout
fi

rm -f constrain.inp xtbtopo* xtbrestart xtbopt* charges *fukui *omega *gfn1 wbo xtblast.xyz NOT_CONVERGED
