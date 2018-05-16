EventTime 1.0
Type=ChangeExtFreq
Label=#1#
Population: Exc1
Receptor: AMPA
FreqExt=200
EndEvent


EventTime 1.0
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.3
GaussSTD=0
EndEvent


EventTime 1.0
Type=ChangeExtFreq
Label=#1#
Population: Exc7
Receptor: AMPA
FreqExt=4000
EndEvent


EventTime 1.0
Type=ChangeMembraneNoise
Label=#1#
Population: Exc7
GaussMean=0.3
GaussSTD=0
EndEvent

%--------------------------------------------------

EventTime 6000.00
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl
FileName:MemPot.dat
Type=MemPot
population:AllPopulation
EndOutputFile

FileName:Freq.dat
Type=FiringRate
FiringRateWinodw=1000
PrintStep=1
population:AllPopulation
EndOutputFile

FileName:Spike.dat
Type=Spike
population:AllPopulation
EndOutputFile

FileName:IOvf.dat
Type=IsynOverLimit
population:AllPopulation
EndOutputFile


EndOutControl
