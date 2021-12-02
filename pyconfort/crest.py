#!/usr/bin/python
from __future__ import print_function, absolute_import

#######################################################################
# Runs crest on an xyz file, can add in options to the code if needed #
#######################################################################

#Python Libraries
import sys, os
import numpy as np
import glob
from optparse import OptionParser
import subprocess
import rdkit
from rdkit.Chem import AllChem

def xyzall_2_xyz(xyzin,name):
	#converting multiple xyz to single
	command_run_1 = ['obabel', xyzin, '-oxyz', '-O'+name+'_conf_.xyz','-m']
	subprocess.run(command_run_1)

def crest_opt(mol, name, dup_data, dup_data_idx, sdwriter, args, log):

	''' run xtb using shell script and args
	Parameters
	----------
	mol : RDKit Mol
		RDkit mol object
	name : str
		Name of molecule
	dup_data : pandas df
		csv
	dup_data_idx : int
		index of pandas df
	sdwriter : Rdkit writer
		sdwriter to write sdf file
	args : pyconfot args
		arguments for crest
	log : Logger
		Logging
	Returns
	-------
	performs crest/cregen conformer search
	'''

	xyzin = name+'.xyz'

	AllChem.EmbedMolecule(mol)
	rdkit.Chem.rdmolfiles.MolToXYZFile(mol, xyzin)

	xyzoutall = name+'_conformers.xyz'
	xyzoutbest = name+'_best.xyz'
	charge = args.charge[0]
	mult = args.mult
	cregen = args.cregen
	cregen_ethr = args.cregen_ethr
	cregen_rthr = args.cregen_rthr
	cregen_bthr = args.cregen_bthr
	cregen_ewin = args.cregen_ewin

	# check whether job has already been run
	if os.path.exists(xyzoutall):
		print('   {0} already exists: skipping crest search'.format(xyzoutall))
		# pass
	else:
		command = [os.path.abspath(os.path.dirname(__file__))+'/run_crest.sh', xyzin, '--xyzoutall', str(xyzoutall),'--xyzoutbest', str(xyzoutbest), '--charge', str(charge)]
		crest_result = subprocess.run(command)

	if cregen:
		xyzcregenensemble = xyzoutall.split('.xyz')[0]+'_cregen.xyz'
		xyzcregenbest = xyzoutbest.split('.xyz')[0]+'_cregen.xyz'

		# check whether job has already been run
		if os.path.exists(xyzcregenensemble):
			print('   {0} already exists: skipping cregen search'.format(xyzcregenensemble))
			# pass
		else:
			command = [os.path.abspath(os.path.dirname(__file__))+'/run_cregen.sh', xyzoutall, xyzoutbest, '--xyzout', str(xyzcregenensemble), '--charge', str(charge),'--ethr', str(cregen_ethr),'--rthr', str(cregen_rthr),'--bthr', str(cregen_bthr),'--ewin', str(cregen_ewin)]
			cregen_result = subprocess.run(command)

			#converting multiple xyz to single
			command_run_1 = ['obabel', xyzcregenensemble, '-oxyz', '-O'+xyzcregenbest,'-l','1']
			subprocess.run(command_run_1)

	if cregen:
		xyzall_2_xyz(xyzcregenensemble,name)
	else:
		xyzall_2_xyz(xyzoutall,name)

	xyz_files = glob.glob(name+'_conf_*.xyz')
	for file in xyz_files:
		name=file.split('.xyz')[0]
		command_xyz = ['obabel', '-ixyz', file, '-osdf', '-O'+name+'.sdf']
		subprocess.call(command_xyz)

	sdf_files = glob.glob(name+'*.sdf')
	for file in sdf_files:
		mol = rdkit.Chem.SDMolSupplier(file)
		sdwriter.write(mol[0])

	dup_data.at[dup_data_idx,'crest-conformers'] = len(xyz_files)

	for f in glob.glob(name+'*_conf_*.xyz'):
		os.remove(f)
	for f in glob.glob(name+'*_conf_*.sdf'):
		os.remove(f)
	for f in glob.glob(name+'*.xyz'):
		os.remove(f)

	return 1
