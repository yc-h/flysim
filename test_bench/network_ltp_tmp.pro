
EventTime 1.0
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.508
GaussSTD=0
EndEvent

EventTime 1.0
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.508
GaussSTD=0
EndEvent

%--------------------------------------------------

EventTime 1800000.0
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl
FileName:MemPot.dat
Type=MemPot
population:AllPopulation
EndOutputFile

FileName:Freq_0.dat
Type=FiringRate
FiringRateWindow=50
PrintStep=10
population:AllPopulation
EndOutputFile

FileName:Spike.dat
Type=Spike
population:AllPopulation
EndOutputFile

FileName:SynWgt.dat
Type=SynapticWeight
population:AllPopulation
EndOutputFile


EndOutControl
