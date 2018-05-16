
EventTime 1.0
Type=ChangeExtFreq
Population: Exc1
Receptor=NMDA
FreqExt=1500
EndEvent


%--------------------------------------------------

EventTime 16000.00
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
