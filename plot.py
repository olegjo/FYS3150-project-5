import matplotlib.pyplot as plt
import numpy as np
import sys

def plotandsave():
	x = []
	u = []

	plt.figure(1)

	for i in range(len(files)):
		filename = files[i]

		x.append([])
		u.append([])
		f = open(filename)
		for line in f:
			numbers = line.split()
			if not float(numbers[1]) == 0:
				x[i].append(float(numbers[0]))
				u[i].append(float(numbers[1]))
		f.close()
		x[i] = np.array(x[i])
		u[i] = np.array(u[i])
		plt.plot(x[i], u[i], label='$\mathrm{%s}$' % (filename[7:-remove_back]))

	u_analytic = np.zeros(len(x[0]))
	for n in range(1, 1000):
		u_analytic += -2/(n*np.pi)*np.sin(n*np.pi*x[0])*np.exp(-n**2*np.pi**2*T)
	u_analytic = u_analytic + 1 - x[0]

	plt.plot(x[0], u_analytic, label='$\mathrm{Analytic}$')
	plt.xlabel('$x \ \ \mathrm{[Rel.\ units]}$', fontsize=18)
	plt.ylabel('$u\ \ \mathrm{[Rel.\ units]}$', fontsize=18)


	plt.legend(loc='best', fontsize=16)

	plt.text(0.2,0.2,'$\Delta x = %s$\n$T=%s$'%(dx,str(T)), fontsize=20,bbox=dict(facecolor='white', alpha=1))

	# finding the errors by least squares
	for method in range(len(files)):
		# err = 0
		# for i in range(len(x)):
		# 	err += (u_analytic[i] - u[method][i])**2

		print files[method], '%.2e' % (np.std(np.abs(u_analytic - u[method]))*100)

	plt.figure(2)

	for method in range(len(files)):
		plt.plot(x[method], abs(np.array(u_analytic) - np.array(u[method]))/(np.array(u_analytic)), label='$\mathrm{%s}$'%str(files[method][7:-remove_back]))

	plt.xlabel('$x\ \ \mathrm{[Rel. \ units]}$', fontsize=18)
	plt.ylabel('$\mathrm{Absolute \ Relative \ error \ \ [\%]}$', fontsize=18)
	plt.legend(loc='best', fontsize = 16)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.ylim([-plt.ylim()[1]*0.05, plt.ylim()[1]])


	plt.figure(1); plt.savefig(filename1)
	plt.close()
	plt.figure(2); plt.savefig(filename2)
	plt.close()

	# plt.show()



files = ['task_d/BackwardEuler_0.1_curved.txt', 'task_d/CrankNickolson_0.1_curved.txt', 'task_d/ForwardEuler_0.1_curved.txt']
remove_back = 15
filename1 = 'FYS3150_project_5_differential_dx01_curved.pdf'
filename2 = 'FYS3150_project_5_differential_error_dx01_curved.pdf'
T = 0.1
dx = '0.1'
plotandsave()

files = ['task_d/BackwardEuler_0.01_curved.txt', 'task_d/CrankNickolson_0.01_curved.txt', 'task_d/ForwardEuler_0.01_curved.txt']
remove_back = 16
filename1 = 'FYS3150_project_5_differential_dx001_curved.pdf'
filename2 = 'FYS3150_project_5_differential_error_dx001_curved.pdf'
T = 0.1
dx = '0.01'
plotandsave()

files = ['task_d/BackwardEuler_0.1_linear.txt', 'task_d/CrankNickolson_0.1_linear.txt', 'task_d/ForwardEuler_0.1_linear.txt']
remove_back = 15
filename1 = 'FYS3150_project_5_differential_dx01_linear.pdf'
filename2 = 'FYS3150_project_5_differential_error_dx01_linear.pdf'
T = 0.3
dt = '0.1'
plotandsave()

files = ['task_d/BackwardEuler_0.01_linear.txt', 'task_d/CrankNickolson_0.01_linear.txt', 'task_d/ForwardEuler_0.01_linear.txt']
remove_back = 16
filename1 = 'FYS3150_project_5_differential_dx001_linear.pdf'
filename2 = 'FYS3150_project_5_differential_error_dx001_linear.pdf'
T = 0.3
dx = '0.01'
plotandsave()














