#!/usr/bin/env python

#####################################################.
# 		   This file stores all the functions 	    #
# 	 used for genrating all nics input/output       #
#####################################################.

from numpy import *
import subprocess, sys, os, math

def find_centroid(ringatoms,CARTESIANS):
		xtot = 0; xvals=[]; yvals=[]; zvals=[]
		for x in ringatoms:
				#print "CARTS", fileData.CARTESIANS[x]
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


def find_coeffplane(ringatoms, CARTESIANS):
		rotated = 0
		xvals, yvals, zvals, xav, yav, zav = find_centroid(ringatoms,CARTESIANS)
		#print xvals, yvals, zvals
		xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum = get_squares_list(ringatoms, xvals, yvals, zvals)
		if xsum == 0.0 and ysum == 0.0:
			rotated = 3
			print("Can't define a ring by points in a line")
			print("This is going to go horribly wrong")
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


def update_coord(NATOMS,ATOMTYPES,CARTESIANS,args,log,name,w_dir_initial):

	#find the ring atoms in the File
	filelines =  open(w_dir_initial+'/'+args.nics_atoms_file,'r').readlines()
	print(filelines)
	ringatoms = []
	for line in (filelines):
		split_line = line.strip().split(',')
		print(split_line)
		print(name)
		print(split_line[0], name)
		if split_line[0] == name.strip():
			for i in range(1,len(split_line)):
				ringatoms.append(int(split_line[i]))
			break

	print(ringatoms)

	coeffplane, xav, yav, zav, rotated = find_coeffplane(ringatoms,CARTESIANS)

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
