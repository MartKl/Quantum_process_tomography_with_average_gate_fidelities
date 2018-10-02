# Quantum_process_tomography_from_average_gate_fidelities
Numerics accompanying the paper Phys. Rev. Lett. [in production, to be updated later], [arXiv:1803.00572](https://arxiv.org/abs/1803.00572). 

We use the shell scripts runmatlab.sh to run the simulations. In this way we can also work with a memory leak of Matlab+cvx+Mosek. 

The simulations are then automatically initialized with init_calc.m and run with calc.m. 
The simulations save the generated data in dataBerlin.mat. 

The data /Haar*/dataBerlin.mat is merged with merge*.m. 
All four different simulations are plotted with plot*.m. 

Four different simulations: Delta(m) and Delta(eta), each for Haar random measurements and random Clifford measurements. 
