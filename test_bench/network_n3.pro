EventTime 1.0
Type=ChangeExtFreq
Label=#1#
Population: Exc1
Receptor: AMPA
FreqExt=650
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
Population: Exc2
Receptor: AMPA
FreqExt=650
EndEvent


EventTime 1.0
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.3
GaussSTD=0
EndEvent




%--------------------------------------------------

EventTime 20000.00
Type=EndTrial
Label=End_of_the_trial
EndEvent

FileName:Freq_0.dat
Type=FiringRate
FiringRateWindow=1000
PrintStep=10
population:AllPopulation
EndOutputFile

FileName:Spike.dat
Type=Spike
population:AllPopulation
EndOutputFile

EndOutControl
