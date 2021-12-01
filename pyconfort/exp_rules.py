import pybel 
from pyconfort.utils import periodic_table

def passes_custom_rules(mol,args,log):
	offset = args.angle_off
	for rule in args.exp_rules:
		passing = apply_rule(rule,mol,log,offset)
		if not passing:
			return False
	return True
def apply_rule(rule,mol,log,offset=0.0):
	if rule == 'Ir_bidentate_x3':
		return apply_Ir_rule(mol,offset)
	return apply_general_rule(rule,mol,log,offset)

def apply_general_rule(rule,mol,log,offset=0.0):
	try:
		syms = rule.split(',')[0].split('-')
		target_angle = float(rule.split(',')[1])
	except IndexError:
		log.write('x  The exp_rules parameter(s) was not correctly defined, this filter will be turned off')
		raise ValueError(f'rule "{rule}" is not properly defined')
	atnums = [periodic_table.index(sym) for sym in syms]
	smarts = pybel.Smarts('[#{}]~[#{}]~[#{}]'.format(*atnums))
	matches = smarts.findall(mol)
	if len(matches) >= 2: 
		log.write(f'x  There are multiple options in exp_rules for {mol.title}, this filter will be turned off')
		log.write(f'x  {mol.title} contain more than one atom that meets the exp_rules criteria, this filter will be turned off')
	
	atoms = [mol.atoms[i].OBAtom for i in matches[0]]
	angle = mol.OBMol.GetAngle(*atoms)
	min_angle = target_angle - offset
	max_angle = target_angle + offset
	if min_angle <= angle <= max_angle:
		return False
	return True
def apply_Ir_rule(mol,offset=0.0):
	smarts = pybel.Smarts('[#77]1~[#7]~*~*~[#6]1')
	matches = smarts.findall(mol)
	min_angle = 180 - offset
	max_angle = 180 + offset
	if len(matches) == 3: 
		# If 3 Py_Ph are coordinated Ns must not be at ~180º
		Ir = mol.atoms[matches[0][0]].OBAtom
		n1 = mol.atoms[matches[0][1]].OBAtom
		n2 = mol.atoms[matches[1][1]].OBAtom
		n3 = mol.atoms[matches[2][1]].OBAtom
		angles_atoms = [(n1,Ir,n2),
						(n1,Ir,n3),
						(n2,Ir,n3)]
		for i,j,k in angles_atoms: 
			angle = mol.OBMol.GetAngle(i,j,k)
			if min_angle <= angle <= max_angle:
				return False

	elif len(matches) == 2:
		# If only 2 Py_Ph are coordinated Ns must be at ~180º
		Ir = mol.atoms[matches[0][0]].OBAtom
		n1 = mol.atoms[matches[0][1]].OBAtom
		n2 = mol.atoms[matches[1][1]].OBAtom
		angle = mol.OBMol.GetAngle(n1,Ir,n2)
		if not(min_angle <= angle <= max_angle):
			return False

	return True