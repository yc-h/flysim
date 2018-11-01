# flysim
This software is a neuron and neural network simulator designed by Yu-Chi in CCLo Lab @2018  
At first, you need GCC 4.8.2 or later version to compile this program.  
Under the same directory of Makefile, please key in:  
  
make  
...  
...  
...  
to compile all files and generate "flysim.out"  


# How to use
flysim -option parameters:  
-pro network.pro     # read protocal file: default=network.pro
-conf network.conf   # read configuration file: default=network.conf
-om my_membrane.dat  # for batch operation: output membrane potential file
-os my_spike.dat     # for batch operation: output spikes file
-or my_rate.dat      # for batch operation: output firing rate: default rate window=50ms, print out=100ms
-rp 1                # set repeat times: default=1
-t 4                 # set multithreading: default=1
-dt 0.1              # time step(default=0.1ms)
-STP                 # use short term plasticity synapse
-STD                 # use short term depression synapse and this option is disabled when -STP used
-LTP                 # use long term plasticity synapse(STDP)
-nmodel GNL2         # neuron model: classical leaky integrate and fire model

  
  
example:  
./flysim.out -conf test_bench/network_n2.conf -pro test_bench/network_n2.pro -nmodel GNL2 -STP  

Configuration file statsitic  
Configuration file setup:100%        
protocal Macros setup  
protocal file setup  
write configuration information to ConfsInfo.log, please wait!  
calculate synapses size:100%        
setup network:100%        
time step:0.1, total neurons=2, total synapses=1  
thread number:1  
write log information to network.log , please wait!  
start processing----------------------------------------------------  
Iteration:1  
Expertment time=15.9999s real time=0.4512544s 100%  


# User guide
File "UserGuide.pdf" is used for simulation

# License
This software is under GPL version 2 or later versions license.
