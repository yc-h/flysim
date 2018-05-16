EventTime 1.0
Type=ChangeExtFreq
Label=#1#
Population: Exc1
Receptor: AMPA
FreqExt=500.0
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
Population: Inh1
Receptor: AMPA
FreqExt=500.0
EndEvent


EventTime 1.0
Type=ChangeMembraneNoise
Label=#1#
Population: Inh1
GaussMean=0.3
GaussSTD=0
EndEvent


%--------------------------------------------------

EventTime 60000.0
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl

FileName:Freq_0.dat
Type=FiringRate
FiringRateWindow=50
PrintStep=10
population:AllPopulation
EndOutputFile

EndOutControl
