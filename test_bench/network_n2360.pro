EventTime 200.0
Type=ChangeExtFreq
Label=#1#
Population: DecisionExcBg
Receptor: AMPA
FreqExt=1300.0
EndEvent

EventTime 200.0
Type=ChangeExtFreq
Label=#1#
Population: Decision1
Receptor: AMPA
FreqExt=1300.0
EndEvent

EventTime 200.0
Type=ChangeExtFreq
Label=#1#
Population: Decision2
Receptor: AMPA
FreqExt=1300.0
EndEvent

EventTime 200.0
Type=ChangeExtFreq
Label=#1#
Population: DecisionInh
Receptor: AMPA
FreqExt=1300.0
EndEvent


%--------------------------------------------------

EventTime 4000.00
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl

FileName:FRates.dat
Type=FiringRate
FiringRateWindow=50
PrintStep=10
population:AllPopulation
EndOutputFile

EndOutControl

