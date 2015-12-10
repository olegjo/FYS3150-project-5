import numpy as np
import matplotlib.pyplot as plt
import os

run_again = False

def histplot(filename, T):
	nZeros = 0
	nOnes = 0
	x = []
	f = open(filename)
	for line in f:
		numbers = line.split()
		if numbers[0] == '#':
			dx = float(numbers[1])
			continue
		if float(numbers[0]) == 0:
			nZeros += 1
			continue
		if float(numbers[0]) > (1-1e-9):
			nOnes += 1
		x.append(float(numbers[0]))
	f.close()
	print '\n', filename,'\nParticles at 0:', nZeros, '\nParticles at 1:', nOnes
	dx = dx*0.2
	X = 0
	bins = []
	while X <= 1:
		bins.append(X)
		X += dx
	plt.figure('dummy')
	data, bins, patches = plt.hist(x, bins)
	if data[0] == 0:
		data[0] = nZeros
	data = data/data[0]
	plt.close('dummy')

	u_analytic = np.zeros(1000)
	x = np.linspace(0, 1, 1000)
	for n in range(1, 1000):
		u_analytic += -2/(n*np.pi)*np.sin(n*np.pi*x)*np.exp(-n**2*np.pi**2*T)

	plt.figure(filename)
	plt.plot(x, u_analytic+1-x, 'r', label='$\mathrm{Analytic}$')
	plt.xlabel('$x \ \ \mathrm{[Rel.\ units]}$', fontsize=18)
	plt.ylabel('$u\ \ \mathrm{[Rel.\ units]}$', fontsize=18)


	plt.bar(bins[:-1], data, width=dx, label='$\mathrm{Numerical \ Random \ walk}$')
	plt.legend(loc='best', fontsize=16)




T = [0.1, 0.3]

if run_again:
	for tStop in T:
		filename = 'RandomWalk_constant_step_%f.txt' % tStop
		os.system('task_f/task_f %f %s' % (tStop, filename))
		filename = 'RandomWalk_gaussian_step_%f.txt' % tStop
		os.system('task_g/task_g %f %s' % (tStop, filename))


for tStop in T:
	filename = 'task_f/RandomWalk_constant_step_%f.txt' % tStop
	histplot(filename, tStop)
	plt.text(0.8, 0.6, '$T=%s$'%str(tStop), fontsize=20,bbox=dict(facecolor='white', alpha=1))
	fig_filename = 'RandomWalk_constant_step_T0%s.pdf' % (str(tStop)[-1])
	plt.savefig(fig_filename)

T = [0.1, 0.3, 0.5]

for tStop in T:
	filename = 'task_g/RandomWalk_gaussian_step_%f.txt' % tStop
	histplot(filename, tStop)
	plt.text(0.8, 0.6, '$T=%s$'%str(tStop), fontsize=20,bbox=dict(facecolor='white', alpha=1))
	fig_filename = 'RandomWalk_gaussian_step_T0%s.pdf' % (str(tStop)[-1])
	plt.savefig(fig_filename)


# plt.show()



















