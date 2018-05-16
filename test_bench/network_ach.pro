DefineMacro
GroupName:Sti2Inh
GroupMembers:Exc6,Exc3
EndGroupMembers

GroupName:Sti3Inh
GroupMembers:Exc5,Exc4,Exc2
EndGroupMembers

EndDefineMacro

EventTime 7.0
Type=ChangeExtFreq
Label=#1#
Population: Sti3Inh
Receptor: AMPA
FreqExt=400
EndEvent


EventTime 3.0
Type=ChangeExtFreq
Label=#1#
Population: Sti2Inh
Receptor: AMPA
FreqExt=500
EndEvent


EventTime 1.0
Type=ChangeExtFreq
Label=#1#
Population: Exc1
Receptor: AMPA
FreqExt=300
EndEvent


EventTime 1.0
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.3
GaussSTD=0
EndEvent

EventTime 11.0
Type=ChangeExtFreq
Label=#1#
Population: AllPopulation
Receptor: AMPA
FreqExt=300
EndEvent



%--------------------------------------------------

EventTime 6000.00
Type=EndTrial
Label=End_of_the_trial
EndEvent

OutControl
FileName:Freq.dat
Type=FiringRate
FiringRateWindow=1000
PrintStep=1
population:AllPopulation
EndOutputFile

FileName:MemPot.dat
Type=MemPot
population:AllPopulation
EndOutputFile

FileName:Spike.dat
Type=Spike
population:AllPopulation
EndOutputFile

FileName:SynGat.dat
Type=GatingVariable
population:AllPopulation
EndOutputFile


EndOutControl
