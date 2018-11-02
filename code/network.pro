DefineMacro
GroupName:Sti2
GroupMembers:Exc2,Exc1
EndGroupMembers
EndDefineMacro

EventTime 3000.0
Type=ChangeExtFreq
Label=#1#
Population: Exc1
Receptor: AMPA
FreqExt=240.0
EndEvent

EventTime 6500.0
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.43
GaussSTD=0.03
EndEvent

EventTime 10100.00
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl
FileName:MemPot.dat
Type=MemPot
population:AllPopulation
EndOutputFile

FileName:Spikes.dat
Type=Spike
population:Sti2
EndOutputFile

FileName:FRates.dat
Type=FiringRate
FiringRateWinodw=50
PrintStep=10
population:Exc2
EndOutputFile

EndOutControl
