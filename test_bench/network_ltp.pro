
EventTime 100
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 101
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.0
GaussSTD=0
EndEvent

EventTime 100
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 101
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.0
GaussSTD=0
EndEvent


EventTime 150
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 151
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.0
GaussSTD=0
EndEvent

EventTime 150
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 151
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.0
GaussSTD=0
EndEvent


EventTime 200
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 201
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.0
GaussSTD=0
EndEvent

EventTime 200
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 201
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.0
GaussSTD=0
EndEvent


EventTime 250
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 251
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.0
GaussSTD=0
EndEvent

EventTime 250
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 251
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.0
GaussSTD=0
EndEvent


EventTime 300
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 301
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.0
GaussSTD=0
EndEvent

EventTime 300
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 301
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.0
GaussSTD=0
EndEvent


EventTime 350
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 351
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.0
GaussSTD=0
EndEvent

EventTime 350
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 351
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.0
GaussSTD=0
EndEvent


EventTime 400
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 401
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.0
GaussSTD=0
EndEvent

EventTime 400
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 401
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.0
GaussSTD=0
EndEvent


EventTime 450
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 451
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.0
GaussSTD=0
EndEvent

EventTime 450
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 451
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.0
GaussSTD=0
EndEvent


EventTime 500
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 501
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.0
GaussSTD=0
EndEvent

EventTime 500
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 501
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.0
GaussSTD=0
EndEvent


EventTime 550
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 551
Type=ChangeMembraneNoise
Label=#1#
Population: Exc1
GaussMean=0.0
GaussSTD=0
EndEvent

EventTime 550
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=15.0
GaussSTD=0
EndEvent

EventTime 551
Type=ChangeMembraneNoise
Label=#1#
Population: Exc2
GaussMean=0.0
GaussSTD=0
EndEvent



%--------------------------------------------------

EventTime 600.0
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
