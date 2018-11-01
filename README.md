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
-pro network.pro     %read protocal file: default=network.pro  
-conf network.conf   %read configuration file: default=network.conf  
-om my_membrane.dat  %for batch opreation, print out membrane potential file  
-os my_spike.dat     %for batch opreation, print out spikes file  
-or my_rate.dat      %for batch opreation, print out firing rate file(rate window=50ms, print time=100ms)  
-rp 1                %set repeat times: default=1  
-t 4                 %set multithreading: default=1  
-s accurate          %for -nmodel sim06, numerical error level of solver:  
                      accurate(RK4), moderate(improved Eular), rough(default, Eular)  
-daemon port         %(experiment function)flysim as daemon for hanitu project: port number  
-dt 0.1              %time step(default=0.1ms)  
-udfsed 1            %user define random seed:0~2^32-1  
-STP                 %use short term plsticity synapse  
-STD                 %use short term depression synapse and this option is disabled when -STP used  
-LTP                 %use long term plsticity synapse(STDP)  
-nmodel LIF          %neuron model:  
                      sim06: capable mode of sim06_10(old version) leaky integrate and fire model  
                      sim07: capable mode of sim07_21(old version) leaky integrate and fire model  
                      LIF(default): classical leaky integrate and fire model  
  
  
example:  
./flysim.out -conf test_bench/network_n2.conf -pro test_bench/network_n2.pro -nmodel LIF -STP  

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
