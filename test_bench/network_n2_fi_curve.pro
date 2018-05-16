EventTime 1.0
Type=ChangeExtFreq
Label=#1#
Population: Exc2
Receptor: AMPA
FreqExt=0.0
EndEvent

EventTime 2.0
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=-0.083
GaussSTD=0.0
EndEvent

%--------------------------------------------------

EventTime 5100.00
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl

FileName:MemPot.dat
Type=MemPot
population:Exc3
population:Exc5
population:Exc2
population:Exc1
EndOutputFile

FileName:Spikes.dat
Type=Spike
population:AllPopulation
EndOutputFile

FileName:FRates.dat
Type=FiringRate
FiringRateWindow=50
PrintStep=50
population:AllPopulation
EndOutputFile

EndOutControl
