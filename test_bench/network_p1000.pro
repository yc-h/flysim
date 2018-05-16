
EventTime 1.0
Type=ChangeExtFreq
Label=#1#
Population: AllPopulation
Receptor: AMPA
FreqExt=3000
EndEvent




%--------------------------------------------------

EventTime 10000
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl

FileName=Spike.dat
Type=Spike
population=AllPopulation
EndOutputFile


EndOutControl
