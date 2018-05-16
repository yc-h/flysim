EventTime 1.0
Type=ChangeExtFreq
Label=#1#
Population: Exc1
Receptor: Ach
FreqExt=350
EndEvent


%--------------------------------------------------

EventTime 1000.00
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl

FileName:Spike.dat
Type=Spike
population:AllPopulation
EndOutputFile

EndOutControl
