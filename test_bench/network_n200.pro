% --------------------------
EventTime 1.

Type=ChangeMembraneNoise
Label=Noise
Population=AllPopulation
GaussMean=0
GaussSTD=0
EndEvent

Type=ChangeExtFreq
Label=begin
Population=AllPopulation
Receptor=AMPA
FreqExt=1.000000
EndEvent

% --------------------------
EventTime 500.

Type=ChangeExtFreq
Label=start_stimulus@Ext8@500@25
Population=Ext8
Receptor=Ach
FreqExt=25.000000
EndEvent

% --------------------------
EventTime 1000.

Type=ChangeExtFreq
Label=stop_stimulus@Ext8@1000
Population=Ext8
Receptor=Ach
FreqExt=0.000000
EndEvent

% --------------------------
EventTime 1500.

Type=ChangeExtFreq
Label=start_response@1500
Population=AllPopulation
Receptor=AMPA
FreqExt=1.000000
EndEvent

% --------------------------
EventTime 2000.

Type=ChangeExtFreq
Label=stop_response@2000
Population=AllPopulation
Receptor=AMPA
FreqExt=1.000000
EndEvent

% --------------------------
EventTime 2100.

Type=EndTrial
Label=End_of_the_trial
EndEvent


