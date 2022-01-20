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


charge=${charge:-0}
nproc=${nproc:-24}
force_const=${force_const:-1.0}
max_cycle=${max_cycle:-100}
json=${json:-true}
xyzout=${xyzout:-xtbopt.xyz}

# input geometry (xyz)
input=$1
file="${input%.*}"
echo -e "-  RUNNING $file WITH xTB \c"

# user-defined variables
echo -e "\$constrain" >> constrain.inp
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        echo -e $1 $2 "\c"
        if [[ "$param" == *dist* ]]
           then
              echo -e "distance: $2" >> constrain.inp
        fi
        if [[ "$param" == *angle* ]]
           then
              echo -e "angle: $2" >> constrain.inp
        fi
        if [[ "$param" == *dihedral* ]]
           then
              echo -e "dihedral: $2" >> constrain.inp
        fi
   fi
  shift
done

# output
outfile="${xyzout%.*}".out
echo -e "-  outfile $outfile \c"
echo -e "-  outxyz $xyzout \c"

# write constraints to file
echo -e "force constant=$force_const" >> constrain.inp
echo -e "\$opt\nmaxcycle=$max_cycle\n\$end" >> constrain.inp

# run xTB if xyz file supplied
if [ -z "$file" ]
then
   echo "NO INPUT!"
else
   xtb $file.xyz --opt --input constrain.inp -c $charge > $outfile && mv xtbopt.xyz $xyzout
fi

rm -f constrain.inp xtbtopo* xtbrestart xtbopt* charges *fukui *omega *gfn1 wbo xtblast.xyz NOT_CONVERGED
