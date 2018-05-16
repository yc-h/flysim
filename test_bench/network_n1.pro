
EventTime 1.0
Type=ChangeMembraneNoise
Population: Exc1
GaussMean=0.85
GaussSTD=0.3
EndEvent




%--------------------------------------------------

EventTime 10000
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl
FileName:MemPot.dat
Type=MemPot
population:AllPopulation
EndOutputFile


FileName=FRate.dat
Type=FiringRate
population=AllPopulation
FiringRateWindow=100
PrintStep=10
EndOutputFile


FileName=Spike.dat
Type=Spike
population=AllPopulation
EndOutputFile


EndOutControl
