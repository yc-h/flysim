# flysim
At first time, you need GCC 4.8.2 or later version to compile this program.  
In the same directory of Makefile, please key in:  
  
make  
...  
...  
to compile all files and generate "flysim.out"  


# How to use
neuron and neural network simulator designed by Yu-Chi in CCLo Lab @2018  
flysim -option parameters:  
-pro network.pro     %read protocal file: default=network.pro  
-conf network.conf   %read configuration file: default=network.conf  
-om my_membrane.dat  %for batch opreation: output membrane potential file  
-os my_spike.dat     %for batch opreation: output spikes file  
-or my_rate.dat      %for batch opreation: output firing rate: default rate window=50ms, print out=100ms  
-rp 1                %set repeat times: default=1  
-t 4                 %set multithreading: default=1  
-s accurate          %for -nmodel GNL, numerical error level of solver:  
                      accurate(RK4), moderate(improved Eular), rough(default, Eular)  
-daemon port         %Flysim as daemon(experiment): port number  
-dt 0.1              %time step(default=0.1ms)  
-udfsed 1            %user define random seed:0~2^32-1  
-STP                 %use short term plsticity synapse  
-STD                 %use short term depression synapse and this option is disabled when -STP used  
-LTP                 %use long term plsticity synapse(STDP)  
-SodCH               %Sodium channel(experiment), only used in LIF  
-nmodel LIF          %neuron model:  
                      sim06: capable mode of sim06_10 leaky integrate and fire model  
                      sim07: capable mode of sim07_21 leaky integrate and fire model  
                      LIF(default): classical leaky integrate and fire model  
                      HH: Hodgkin-Huxley model  
  
  
examples:  
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


you can get simulation result from flysim

# User guide
File "UserGuide.pdf" write in Chinese for simulator perform simulation

# License
This software is under GPL version 2 or later versions license.
