# flysim
neuron and neural network simulator designed by Yu-Chi in CCLo Lab @2018
flysim -option parameters
examples:
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

# User guide
File "UserGuide.pdf" write in Chinese for simulator perform simulation

# License
This software is under GPL version 2 or later versions license.
