EventTime 1.0
Type=ChangeExtFreq
Label=#1#
Population: AllPopulation
Receptor: AMPA
FreqExt=300
EndEvent

EventTime 1.0
Type=ChangeMembraneNoise
Label=#1#
Population: AllPopulation
GaussMean=0.3
GaussSTD=0
EndEvent



%--------------------------------------------------

EventTime 6000.00
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl

FileName:Spike.dat
Type=Spike
population:AllPopulation
EndOutputFile

EndOutControl
