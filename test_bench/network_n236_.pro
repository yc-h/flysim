
DefineMacro
GroupName:group1
GroupMembers:Decision1,DecisionExcBg,PFC1
EndGroupMembers
EndDefineMacro

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

EventTime 300.0
Type=ChangeExtFreq
Label=#1#
Population: group1
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

EventTime 400.0
Type=ChangeExtFreq
Label=#1#
Population: group1
Receptor: AMPA
FreqExt=1300.0
EndEvent

EventTime 500.0
Type=ChangeExtFreq
Label=#1#
Population: DecisionInh
Receptor: AMPA
FreqExt=1300.0
EndEvent



EventTime 3000.0
Type=Resume
InFile=save.dat
SubType: ResumeAll
Population: DecisionExcBg
EndEvent

EventTime 3000.0
Type=Resume
InFile=save.dat
SubType: ResumeAll
Population: Decision1
EndEvent

EventTime 3000.0
Type=Resume
InFile=save.dat
SubType: ResumeAll
Population: Decision2
EndEvent

EventTime 3000.0
Type=Resume
InFile=save.dat
SubType: ResumeAll
Population: Fixation
EndEvent

EventTime 3000.0
Type=Resume
InFile=save.dat
SubType: ResumeAll
Population: DecisionInh
EndEvent

EventTime 3000.0
Type=Resume
InFile=save.dat
SubType: ResumeAll
Population: PFC1
EndEvent


%--------------------------------------------------

EventTime 4000.00
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl

FileName:Spikes.dat
Type=Spike
population:AllPopulation
EndOutputFile

FileName:MemPot.dat
Type=MemPot
population:AllPopulation
EndOutputFile

FileName:FRates.dat
Type=FiringRate
FiringRateWindow=50
PrintStep=10
population:AllPopulation
EndOutputFile

FileName:Isyn.dat
Type=IsynSeparatedSum
population:group1
EndOutputFile

FileName:save.dat
Type=SaveAllStates
population:AllPopulation
EventTime=1.0
EndOutputFile


FileName:GatVar.dat
Type=GatingVariable
population:PFC1
EndOutputFile

EndOutControl

