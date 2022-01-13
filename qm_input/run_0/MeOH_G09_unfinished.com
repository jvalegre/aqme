%nprocshared=12
%mem=12GB
# opt=calcfc freq=noraman cc-pvtz scrf=(solvent=chloroform,pcm) pbe1pbe g09defaults

MeOH_G09_unfinished

0 1
 C   0.67496200  -0.01912000   0.00000000
 H   1.13330900   0.94774000  -0.00000300
 H   0.97796200  -0.55747600  -0.87365000
 H   0.97796200  -0.55747000   0.87365300
 O  -0.74748200   0.12768700   0.00000000
 H  -1.15914600  -0.73957000   0.00000000

