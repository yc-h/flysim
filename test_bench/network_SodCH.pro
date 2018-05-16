EventTime 1.0
Type=ChangeExtFreq
Label=#1#
Population: Exc1
Receptor: AMPA
FreqExt=0
EndEvent


EventTime 1.0
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=-1.0
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

FileName:Spike.dat
Type=Spike
population:AllPopulation
EndOutputFile

FileName:Freq_0.dat
Type=FiringRate
FiringRateWindow=50
PrintStep=10
population:AllPopulation
EndOutputFile

EndOutControl
