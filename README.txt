calcSpaceStation.m: calculate steady-state results for the analytical Space Station model

calcStopRestart.m: calculate steady-state results for the analytical CoG model

calcSpaceStation_Branch.m: space station model for symmetric branched geometry

calcSpaceStation_antero.m: space station model with only anterograde mitochondria allowed to fuse

calcStopRestart_localtrans.m: CoG model with local translation included


--------------------------------------------
Simulation Code is provided (in fortran) in the simulations subdirectory
--------------------------------------------

Sequence of steps to run and visualize simulations:

1. Compile the code: open a terminal at the simulation code directory (either simulations/sngphago_kp for the CoG model or simulations/ssphago for the SS model)
and compile by running 'make'. The simulation executable is now at simulations/*.exe

2. Specify the parameter file (param.<paramfilename>): See examples/param.example_ss and examples/param.example_cog for some example parameters files to run.

3. Navigate to the directory with your parameter file. Run the compiled code with your chosen parameters as
../sngkpphago.exe <paramfilename> 

For example: 
../sngkpphago.exe example_cog (for the CoG model)
or
../ssphago.exe example_ss (for the SS model)

The simulation snapshots are stored in your parameter directory as <paramfilename>.snap.out 

4. To visualize the snapshots, use the matlab code processSims/process_examples.m
