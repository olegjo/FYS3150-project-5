# FYS3150-project-5
This is the fifth and final numerical project in FYS3150 - Computational Physics at University of Oslo, fall 2015.

## Abstract
In this report, I investigate the properties of the forward Euler, backward Euler and 
Crank-Nicolson schemes and use these methods to solve the 1+1D diffusion equation in 
light of diffusion of neurotransmitters in the synaptic cleft between brain cells. 
I also investigate the method’s stability and find that all produce results with good 
agreement to the analytical solution, with the poorest result being 0.4% off.

I also look at the diffusion process as a series of random walkers walking across 
the synaptic cleft and find that the best results are obtained when the particles 
are constrained to move with a constant step length, as compared to one with a 
normal distribution. In light of these results I propose an explanation of the underlying 
physics of this particular diffusion process and how the density of the diffusive medium 
may alter the results in a way that renders the normally distributed step lengths a more 
feasible approach.

## Dependencies
C++11 with standard vector class

Python 2.7 with numpy and matplotlib

## Overview
The main directory has two important files: ```plot.py``` which produces the plots for the differential methods(FE, BE and CN)
and ```histplot.py``` which produces the histograms for the random walk method. The latter also includes a test for the
random walk methods.

The subdirectories are:
###### task_d
Implementation the differential methods
###### task_f
Implementation of random walk with constant step length
###### task_g
Implementation of random walk with normally distributed step lengths

## Benchmarks
All files to produce the plots are stored, however, if desireable to reproduce the results, the necessary runs are:

```task_d/task_d 0.1 0.1 0.1_curved```

```task_d/task_d 0.1 0.3 0.1_linear```

```task_d/task_d 0.01 0.1 0.01_curved```

```task_d/task_d 0.01 0.3 0.01_linear```


```task_f/task_f 0.1 RandomWalk_constant_step_0.100000.txt```

```task_f/taks_f 0.3 RandomWalk_constant_step_0.300000.txt```

```task_g/task_g 0.1 RandomWalk_gaussian_step_0.100000.txt```

```task_g/task_g 0.3 RandomWalk_gaussian_step_0.300000.txt```

```task_g/task_g 0.5 RandomWalk_gaussian_step_0.500000.txt```




