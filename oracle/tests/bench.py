#creates a folder named ddmyyyy-hh-m-sec and chucks everything in it. The plot is in the same folder, with the same name

import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import glob
import datetime

solver='classical'
yaml='hermes_classical_2.yaml'

names='18Sco,61CygA,61CygB,alfCet,alfTau,Arcturus,betVir,delEri,epsVir,etaBoo,gamSge,Gmb1830,HD107328,HD122563,HD140283,\
HD220009,HD22879,HD84937,muCas,muLeo,Procyon,tauCet'
names=names.split(',')

#sorted by names 
#[(teff,dteff),(logg,dlogg),(feh),(xi)]
bench=[[(5747,39),(4.43,0.01),(0.03),(2.2)],
[(4339,27),(4.43,0.16),(-0.33),(0.0)],
[(4045,25),(4.53,0.04),(-0.38),(1.7)],
[(3796,65),(0.91,0.08),(-0.45),(3)],
[(3927,40),(1.22,0.1),(-0.37),(5)],
[(4247,37),(1.59,0.04),(-0.52),(3.8)],
[(6083,41),(4.08,0.01),(0.24),(2)],
[(5045,65),(3.77,0.02),(0.06),(0.7)],
[(4983,61),(2.77,0.01),(0.15),(2)],
[(6105,28),(3.8,0.02),(0.32),(12.7)],
[(3807,49),(1.05,0.1),(-0.17),(6)],
[(4827,55),(4.6,0.03),(-1.46),(0.5)],
[(4590,59),(2.2,0.07),(-0.33),(1.9)],
[(4608,60),(1.61,0.07),(-2.64),(5)],
[(5720,120),(3.67,0.04),(-2.36),(5)],
[(4266,54),(1.43,0.1),(-0.74),(1)],
[(5789,89),(4.23,0.03),(-0.86),(4.4)],
[(6275,97),(4.11,0.06),(-2.03),(5.2)],
[(5308,29),(4.41,0.02),(-0.81),(0)],
[(4433,60),(2.5,0.07),(0.25),(5.1)],
[(6545,84),(3.99,0.02),(0.01),(2.8)],
[(5331,43),(4.44,0.02),(-0.49),(1.1)]]

date=datetime.datetime.now()
filename=str(date.day)+str(date.month)+str(date.year)+'-'+':'.join([str(date.hour),str(date.minute),str(date.second)])
os.system('mkdir %s' %filename)

os.chdir(filename)
for i in names:
	os.system('oracle solve %s -o %s ../%s ../data/benchmarks/%s/*noresample*.txt'\
	%(solver,i,yaml,i))

logs=[i[:-4] for i in glob.glob('*log')]


o_val=[] #sorted by logs
o_val_names=[]
for i in logs:
	proc = subprocess.Popen("grep 'Optimised parameters' %s.log | tail -1 " \
   	%i, stdout=subprocess.PIPE, shell=True)
   	tmp = filter(None,proc.stdout.read().split(' '))
	try:
		o_val.append([float(tmp[9]),float(tmp[13]),float(tmp[16][:-1]),float(tmp[-2])])
		o_val_names.append(i)
	except ValueError:
		print i
		continue
bench_r=[] #bench reordered to logs

for i in range(len(o_val)):
	bench_r.append(bench[names.index(o_val_names[i])]) 

plt.subplot(2, 2, 1)
plt.errorbar(range(len(bench_r)),[i[0][0] for i in bench_r],yerr= [i[0][1] for i in bench_r],fmt=None,label='Jofre')
plt.scatter(range(len(bench_r)),[i[0] for i in o_val],marker='x',c='r',label='oracle')
plt.xticks(range(len(bench_r)),o_val_names,rotation='vertical')
plt.ylabel('teff(K)')

plt.subplot(2, 2, 2)
plt.errorbar(range(len(bench_r)),[i[1][0] for i in bench_r],yerr= [i[1][1] for i in bench_r],fmt=None)
plt.scatter(range(len(bench_r)),[i[1] for i in o_val],marker='x',c='r')
plt.xticks(range(len(bench_r)),o_val_names,rotation='vertical')
plt.ylabel('logg')

plt.subplot(2, 2, 3)
plt.scatter(range(len(bench_r)),[i[2] for i in bench_r],label='Jofre')
plt.scatter(range(len(bench_r)),[i[2] for i in o_val],marker='x',c='r',label='oracle')
plt.xticks(range(len(bench_r)),o_val_names,rotation='vertical')
plt.ylabel('feh(dex)')

plt.subplot(2, 2, 4)
plt.scatter(range(len(bench_r)),[i[-1] for i in bench_r])
plt.scatter(range(len(bench_r)),[i[-1] for i in o_val],marker='x',c='r')
plt.xticks(range(len(bench_r)),o_val_names,rotation='vertical')
plt.ylabel('xi(km/s)')
plt.savefig(filename)
