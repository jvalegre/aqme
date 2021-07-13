import shutil
import os
from rdkit import Chem

#  Functions brought from qprep

# WRITE SDF FILES FOR xTB AND ANI1
def write_confs(conformers, energies,selectedcids, name, args, program,log):
	if len(conformers) > 0:
		# name = name.split('_'+args.CSEARCH)[0]# a bit hacky
		sdwriter = Chem.SDWriter(name+'_'+program+args.output)

		write_confs = 0
		for cid in selectedcids:
			sdwriter.write(conformers[cid])
			write_confs += 1

		if args.verbose:
			log.write("o  Writing "+str(write_confs)+ " conformers to file " + name+'_'+program+args.output)
		sdwriter.close()
	else:
		log.write("x  No conformers found!")

# MOVES SDF FILES TO THEIR CORRESPONDING FOLDERS
def moving_files(destination,src,file):
	try:
		os.makedirs(destination)
		shutil.move(os.path.join(src, file), os.path.join(destination, file))
	except OSError:
		if  os.path.isdir(destination):
			shutil.move(os.path.join(src, file), os.path.join(destination, file))
		else:
			raise