#!/usr/bin/env python

#####################################################.
#            This file stores all the functions         #
#      used for genrating all nics input/output       #
#####################################################.

from numpy import *
import subprocess, sys, os, math
from pyconfort.qprep_gaussian import moving_files
import pandas as pd
from pyconfort.argument_parser import possible_atoms

possible_atoms = possible_atoms()

def get_coords_normal(outlines, stand_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS):
	for i in range(stand_or+5,stand_or+5+NATOMS):
		massno = int(outlines[i].split()[1])
		if massno < len(possible_atoms):
			atom_symbol = possible_atoms[massno]
		else:
			atom_symbol = "XX"
		ATOMTYPES.append(atom_symbol)
		CARTESIANS.append([float(outlines[i].split()[3]), float(outlines[i].split()[4]), float(outlines[i].split()[5])])
	return ATOMTYPES, CARTESIANS,stand_or


def find_centroid(ringatoms,CARTESIANS):
	xtot = 0; xvals=[]; yvals=[]; zvals=[]
	for x in ringatoms:
		xtot = xtot + CARTESIANS[x][0]
		xvals.append(CARTESIANS[x][0])
	xav = xtot/len(ringatoms)
	ytot = 0
	for x in ringatoms:
		ytot = ytot + CARTESIANS[x][1]
		yvals.append(CARTESIANS[x][1])
	yav = ytot/len(ringatoms)
	ztot = 0
	for x in ringatoms:
		ztot = ztot + CARTESIANS[x][2]
		zvals.append(CARTESIANS[x][2])
	zav = ztot/len(ringatoms)

	#print "Centroid at:", xav, yav, zav  #gives position of centroid
	return xvals, yvals, zvals, xav, yav, zav

def get_squares_list(ringatoms, xvals, yvals, zvals):
####################Necessary summations
	xysum = 0; y2sum = 0; x2sum = 0; zsum = 0; ysum = 0; xsum = 0; xzsum = 0; yzsum = 0
	for n in range(len(ringatoms)):
		xy = xvals[n]*yvals[n]
		xysum = xy+xysum
		xz = xvals[n]*zvals[n]
		xzsum = xz+xzsum
		yz = yvals[n]*zvals[n]
		yzsum = yz+yzsum
		x = xvals[n]
		xsum = x+xsum
		y = yvals[n]
		ysum = y+ysum
		z = zvals[n]
		zsum = z+zsum
		x2 = xvals[n]*xvals[n]
		x2sum = x2+x2sum
		y2 = yvals[n]*yvals[n]
		y2sum = y2+y2sum
	return xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum

def do_matrix_stuff(xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum, ringatoms):
	###################Matrix and vector used for least squares best fit plane
	a=matrix([[x2sum, xysum, xsum],[xysum, y2sum, ysum],[xsum, ysum, len(ringatoms)]]) #3x3 matrix
	b=matrix([[xzsum],[yzsum],[zsum]]) #3x1 matrix
	try: coeffplane=a.I*b
	except linalg.linalg.LinAlgError: coeffplane = matrix([[0.0],[0.0],[0.0]])
	return coeffplane

def find_coeffplane(ringatoms, CARTESIANS,log):
	rotated = 0
	xvals, yvals, zvals, xav, yav, zav = find_centroid(ringatoms,CARTESIANS)
	#print xvals, yvals, zvals
	xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum = get_squares_list(ringatoms, xvals, yvals, zvals)
	if xsum == 0.0 and ysum == 0.0:
		rotated = 3
		log.write("x  Can't define a ring by points in a line")

	if xsum == 0.0:
		new_xvals = yvals
		new_yvals = zvals
		new_zvals = xvals
		xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum = get_squares_list(ringatoms, xvals, yvals, zvals)
		rotated = 1
	if ysum == 0.0:
		new_xvals = zvals
		new_yvals = xvals
		new_zvals = yvals
		xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum = get_squares_list(ringatoms, xvals, yvals, zvals)
		rotated = 2

	coeffplane = do_matrix_stuff(xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum, ringatoms)
	return coeffplane, xav, yav, zav, rotated

def update_coord(NATOMS,ATOMTYPES,CARTESIANS,args,log,name,w_dir_initial,type):
	#find the ring atoms in the File

	filelines =  open(w_dir_initial+'/'+args.nics_atoms_file,'r').readlines()
	ringatoms = []
	for line in (filelines):
		split_line = line.strip().split(',')
		if split_line[0] == name.strip().split()[0]:
			for i in range(1,len(split_line)):
				ringatoms.append(int(split_line[i])-1)
			break

	coeffplane, xav, yav, zav, rotated = find_coeffplane(ringatoms,CARTESIANS,log)

	xcoeff= coeffplane.tolist()[0][0]
	ycoeff= coeffplane.tolist()[1][0]
	cval= coeffplane.tolist()[2][0]

	rawvector=array([xcoeff,ycoeff,-1]) #Need to make into unit vector

	x=float(rawvector[0])
	y=float(rawvector[1])
	z=float(rawvector[2])
	#print x,y,z
	normfactor=1/(x**2+y**2+z**2)**0.5
	x=x*normfactor; y=y*normfactor; z=z*normfactor
	if z<0: z=-z;y=-y;x=-x #Sign flip if z is negative
	#print "Unit vector:", x, y, z #The length of this vector is 1
	#print "NICS 1 point", x+xav, y+yav, z+zav
	if rotated == 1:
		log.write("************ coordinated system was rotated! ***********")
		old_x = z
		old_y = x
		old_z = y
		if old_z<0: old_z=-old_z;old_y=-old_y;old_x=-old_x
		log.write("Unit vector:", old_x, old_y, old_z)
		x = old_x
		y = old_y
		z = old_z
	if rotated == 2:
		log.write("************ coordinated system was rotated! ***********")
		old_x = y
		old_y = z
		old_z = x
		if old_z<0: old_z=-old_z;old_y=-old_y;old_x=-old_x
		log.write("Unit vector:", old_x, old_y, old_z)
		x = old_x
		y = old_y
		z = old_z
	if rotated == 3:
		log.write("didn't I tell you this was a bad idea?")

	if type=='write':
		spacing = float(args.nics_range)/float(args.nics_number)
		for w in range(-args.nics_number,args.nics_number+1):
			scalefactor = w*spacing
			xscale=x*scalefactor
			yscale=y*scalefactor
			zscale=z*scalefactor
			NATOMS += 1
			ATOMTYPES.append("Bq")
			CARTESIANS.append([xav+xscale, yav+yscale, zav+zscale])

		return NATOMS,ATOMTYPES,CARTESIANS
	if type=='read':
		return x,y,z,xav, yav, zav

def getSHIELDING(outlines,args,log):
	NATOMS,stand_or = 0,0
	ATOMTYPES, CARTESIANS,NMR  = [],[],[]
	#find NATOMS
	for i in range(0,len(outlines)):
		if outlines[i].find("Input orientation") > -1:
			stand_or = i
		if outlines[i].find("Rotational constants") > -1 and outlines[i-1].find("-------") > -1:
			NATOMS = i-stand_or-6
			break

	ATOMTYPES, CARTESIANS,stand_or = get_coords_normal(outlines, stand_or, NATOMS, possible_atoms, ATOMTYPES, CARTESIANS)

	NMR = []
	for i in range(0,len(outlines)):
		if outlines[i].find(" SCF GIAO Magnetic shielding tensor (ppm):") > -1:
			j = 0
			while j < NATOMS*5:
				item = {}
				item['atom_index'] = int(outlines[i+1+j].split()[0])
				item['elementID'] = outlines[i+1+j].split()[1]
				item['isotropic'] = float(outlines[i+1+j].split()[4])
				item['anisotropy'] = float(outlines[i+1+j].split()[7])
				item['xx'] = float(outlines[i+2+j].split()[1])
				item['yx'] = float(outlines[i+2+j].split()[3])
				item['zx'] = float(outlines[i+2+j].split()[5])
				item['xy'] = float(outlines[i+3+j].split()[1])
				item['yy'] = float(outlines[i+3+j].split()[3])
				item['zy'] = float(outlines[i+3+j].split()[5])
				item['xz'] = float(outlines[i+4+j].split()[1])
				item['yz'] = float(outlines[i+4+j].split()[3])
				item['zz'] = float(outlines[i+4+j].split()[5])
				item['eigenvalues'] = [float(outlines[i+5+j].split()[1]), float(outlines[i+5+j].split()[2]), float(outlines[i+5+j].split()[3])]
				NMR.append(item)
				j += 5

	return NMR, CARTESIANS,NATOMS,ATOMTYPES

def calculate_nics_parameters(log_files,args,log,w_dir_initial,name_mol,lot,bs):

	total_data = pd.DataFrame()
	dist_data = pd.DataFrame()

	for counter,log in enumerate(log_files):
		name = log.split('.log')[0]

		outlines =  open(log,'r').readlines()
		NMR, CARTESIANS,NATOMS,ATOMTYPES = getSHIELDING(outlines,args,log)
		x,y,z,xav, yav, zav = update_coord(NATOMS,ATOMTYPES,CARTESIANS,args,log,name_mol,w_dir_initial,'read')

		total_data.at[counter,'Name'] = name_mol
		total_data.at[counter,'log'] = name
		dist_data.at[counter,'Name'] = name_mol
		dist_data.at[counter,'log'] = name

		for tensor in NMR:
			if tensor['elementID'] == "Bq":
				xdist = CARTESIANS[(tensor['atom_index']-1)][0] - xav
				ydist = CARTESIANS[(tensor['atom_index']-1)][1] - yav
				zdist = CARTESIANS[(tensor['atom_index']-1)][2] - zav
				dist = xdist*x + ydist*y + zdist*z

				nics_iso = tensor['isotropic']*-1
				nics_oop = (tensor['xx']*x+tensor['yy']*y+tensor['zz']*z)*-1
				nics_ip = (nics_iso * 3 - nics_oop)/2

				total_data.at[counter,'iso-Bq-'+str(tensor['atom_index'])] = nics_iso
				total_data.at[counter,'oop-Bq-'+str(tensor['atom_index'])] = nics_oop
				total_data.at[counter,'ip-Bq-'+str(tensor['atom_index'])] = nics_ip

				dist_data.at[counter,str(tensor['atom_index'])] = dist


	#creating folder for all molecules to write geom parameter
	if str(bs).find('/') > -1:
		folder = w_dir_initial + '/QPRED/nics/all_confs_nics/'+str(lot)+'-'+str(bs).split('/')[0]
	else:
		folder = w_dir_initial + '/QPRED/nics/all_confs_nics/'+str(lot)+'-'+str(bs)
	try:
		os.makedirs(folder)
	except OSError:
		if os.path.isdir(folder):
			pass

	total_data.to_csv(folder+'/'+name_mol+'-all-nics-data.csv',index=False)
	dist_data.to_csv(folder+'/'+name_mol+'-all-nics-dist-data.csv',index=False)

def calculate_boltz_for_nics(val,args,log,name,w_dir_fin,w_dir_initial,lot,bs):

	# GoodVibes must be installed as a module (through pip or conda)
	cmd_boltz = ['python','-m', 'goodvibes', '--boltz', '--output', name ]
	for file in val:
		cmd_boltz.append(file)
	subprocess.call(cmd_boltz)

	#writing to coorect places
	if str(bs).find('/') > -1:
		destination = w_dir_initial+'/QPRED/nics/boltz/'+str(lot)+'-'+str(bs).split('/')[0]
	else:
		destination = w_dir_initial+'/QPRED/nics/boltz/'+str(lot)+'-'+str(bs)
	moving_files(destination, os.getcwd(),'Goodvibes_'+name+'.dat')


def calculate_avg_nics(val,args,log,name,w_dir_fin,w_dir_initial,lot,bs):
	if str(bs).find('/') > -1:
		nics_file = w_dir_initial + '/QPRED/nics/all_confs_nics/'+str(lot)+'-'+str(bs).split('/')[0]+'/'+name+'-all-nics-data.csv'
		dist_file = w_dir_initial + '/QPRED/nics/all_confs_nics/'+str(lot)+'-'+str(bs).split('/')[0]+'/'+name+'-all-nics-dist-data.csv'
		file = w_dir_initial+'/QPRED/nics/boltz/'+str(lot)+'-'+str(bs).split('/')[0]+'/Goodvibes_'+name+'.dat'
	else:
		nics_file = w_dir_initial + '/QPRED/nics/all_confs_nics/'+str(lot)+'-'+str(bs)+'/'+name+'-all-nics-data.csv'
		dist_file = w_dir_initial + '/QPRED/nics/all_confs_nics/'+str(lot)+'-'+str(bs)+'/'+name+'-all-nics-dist-data.csv'
		file = w_dir_initial+'/QPRED/nics/boltz/'+str(lot)+'-'+str(bs)+'/Goodvibes_'+name+'.dat'

	df_nics =  pd.read_csv(nics_file)
	df_dist =  pd.read_csv(dist_file)
	outlines = open(file,"r").readlines()

	#reading the data from boltz file
	for i in range(len(outlines)):
		# I remove the NMR from the file names using [0:-4]
		if outlines[i].find('   ***************************************************************************************************************************************\n') > -1 and outlines[i-1].find('   Structure') > -1:
			start_line = i+1
		elif outlines[i].find('   ***************************************************************************************************************************************\n') > -1:
			end_line = i

	boltz_values = pd.DataFrame()
	counter = 0
	for i in range(start_line,end_line):
		values = outlines[i].split()
		boltz_values.at[counter,'log'] = values[1]+'_nics'
		boltz_values.at[counter,'boltz'] = values[-1]
		counter+=1

	# multiply actual file and sum and add write to new FIL
	df_nics = df_nics.drop(['Name'], axis=1)
	df_dist = df_dist.drop(['Name'], axis=1)

	df_nics = df_nics.set_index('log')
	df_dist = df_dist.set_index('log')
	boltz_values =  boltz_values.set_index('log')

	avg_data_iso = pd.DataFrame()
	avg_data_oop = pd.DataFrame()
	avg_data_ip = pd.DataFrame()
	for col in df_nics.columns:
		df_nics[col] = df_nics[col].astype(float)*(boltz_values['boltz'].astype(float))
		if col.split('-')[0]=='oop':
			avg_data_oop.at[col,'oop'] = df_nics[col].sum()
			for col2 in df_dist.columns:
				if col.split('-')[2]==col2:
					avg_data_oop.at[col,'dist'] = df_dist[col2].mean()
		if col.split('-')[0]=='ip':
			avg_data_ip.at[col,'ip'] = df_nics[col].sum()
			for col2 in df_dist.columns:
				if col.split('-')[2]==col2:
					avg_data_ip.at[col,'dist'] = df_dist[col2].mean()
		if col.split('-')[0]=='iso':
			avg_data_iso.at[col,'iso'] = df_nics[col].sum()
			for col2 in df_dist.columns:
				if col.split('-')[2]==col2:
					avg_data_iso.at[col,'dist'] = df_dist[col2].mean()

	avg_data_oop['Atom-Number'] = avg_data_oop.index
	avg_data_ip['Atom-Number'] = avg_data_ip.index
	avg_data_iso['Atom-Number'] = avg_data_iso.index
	if str(bs).find('/') > -1:
		folder = w_dir_initial + '/QPRED/nics/average_nics/'+str(lot)+'-'+str(bs).split('/')[0]
	else:
		folder = w_dir_initial + '/QPRED/nics/average_nics/'+str(lot)+'-'+str(bs)
	try:
		os.makedirs(folder)
	except OSError:
		if os.path.isdir(folder):
			pass

	avg_data_oop.to_csv(folder+'/'+name+'-avg-nics-oop.csv',index=False)
	avg_data_ip.to_csv(folder+'/'+name+'-avg-nics-ip.csv',index=False)
	avg_data_iso.to_csv(folder+'/'+name+'-avg-nics-iso.csv',index=False)
