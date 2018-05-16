%@protocol@0-0
%@input@[0, 100]
%@ring@[0, 100]
%@input_region@[4]
%@start_stimulus@100
%@stop_stimulus@1100
%@end_trial@5000

DefineMacro
GroupName:Ring_Neuron
%group=34
GroupMembers:34
EndGroupMembers

GroupName:_macro_0
%group=_region_PB-L4 _region_PB-R5
GroupMembers:21,29
EndGroupMembers

GroupName:_region_EB-L1C
%group=4 5
GroupMembers:4,5
EndGroupMembers

GroupName:_region_EB-L2C
%group=13 14
GroupMembers:13,14
EndGroupMembers

GroupName:_region_EB-L3C
%group=5 6
GroupMembers:5,6
EndGroupMembers

GroupName:_region_EB-L4C
%group=14 15
GroupMembers:14,15
EndGroupMembers

GroupName:_region_EB-L5C
%group=6 7
GroupMembers:6,7
EndGroupMembers

GroupName:_region_EB-L6C
%group=15 16
GroupMembers:15,16
EndGroupMembers

GroupName:_region_EB-L7C
%group=7 8
GroupMembers:7,8
EndGroupMembers

GroupName:_region_EB-L8C
%group=9 16 17
GroupMembers:9,16,17
EndGroupMembers

GroupName:_region_EB-R1C
%group=12 13
GroupMembers:12,13
EndGroupMembers

GroupName:_region_EB-R2C
%group=3 4
GroupMembers:3,4
EndGroupMembers

GroupName:_region_EB-R3C
%group=11 12
GroupMembers:11,12
EndGroupMembers

GroupName:_region_EB-R4C
%group=2 3
GroupMembers:2,3
EndGroupMembers

GroupName:_region_EB-R5C
%group=10 11
GroupMembers:10,11
EndGroupMembers

GroupName:_region_EB-R6C
%group=1 2
GroupMembers:1,2
EndGroupMembers

GroupName:_region_EB-R7C
%group=9 10
GroupMembers:9,10
EndGroupMembers

GroupName:_region_EB-R8C
%group=0 1 8
GroupMembers:0,1,8
EndGroupMembers

GroupName:_region_EB_C
%group=10 7 15 14 1 8 13 5 12 0 2 11 4 9 17 6 3 16
GroupMembers:0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17
EndGroupMembers

GroupName:_region_EB_P
%group=
GroupMembers:
EndGroupMembers

GroupName:_region_PB
%group=22 25 20 19 32 26 31 21 18 29 28 27 33 24 30 23
GroupMembers:18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33
EndGroupMembers

GroupName:_region_PB-L1
%group=26
GroupMembers:26
EndGroupMembers

GroupName:_region_PB-L2
%group=27
GroupMembers:27
EndGroupMembers

GroupName:_region_PB-L3
%group=28
GroupMembers:28
EndGroupMembers

GroupName:_region_PB-L4
%group=29
GroupMembers:29
EndGroupMembers

GroupName:_region_PB-L5
%group=30
GroupMembers:30
EndGroupMembers

GroupName:_region_PB-L6
%group=31
GroupMembers:31
EndGroupMembers

GroupName:_region_PB-L7
%group=32
GroupMembers:32
EndGroupMembers

GroupName:_region_PB-L8
%group=33
GroupMembers:33
EndGroupMembers

GroupName:_region_PB-R1
%group=25
GroupMembers:25
EndGroupMembers

GroupName:_region_PB-R2
%group=24
GroupMembers:24
EndGroupMembers

GroupName:_region_PB-R3
%group=23
GroupMembers:23
EndGroupMembers

GroupName:_region_PB-R4
%group=22
GroupMembers:22
EndGroupMembers

GroupName:_region_PB-R5
%group=21
GroupMembers:21
EndGroupMembers

GroupName:_region_PB-R6
%group=20
GroupMembers:20
EndGroupMembers

GroupName:_region_PB-R7
%group=19
GroupMembers:19
EndGroupMembers

GroupName:_region_PB-R8
%group=18
GroupMembers:18
EndGroupMembers

GroupName:_type_EIP
%group=0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
GroupMembers:0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17
EndGroupMembers

GroupName:_type_PEI
%group=18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33
GroupMembers:18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33
EndGroupMembers

EndDefineMacro

% --------------------------
EventTime 0.
Type=ChangeExtFreq
Label=setup
Population=AllPopulation
Receptor=Ach
FreqExt=0
EndEvent

% --------------------------
EventTime 100.
Type=ChangeExtFreq
Label=turn_on
Population=Ring_Neuron
Receptor=Ach
FreqExt=100
EndEvent

EventTime 100.
Type=ChangeExtFreq
Label=turn_on
Population=_macro_0
Receptor=Ach
FreqExt=100
EndEvent

% --------------------------
EventTime 1100.
Type=ChangeExtFreq
Label=turn_off
Population=Ring_Neuron
Receptor=Ach
FreqExt=0
EndEvent

EventTime 1100.
Type=ChangeExtFreq
Label=turn_off
Population=_macro_0
Receptor=Ach
FreqExt=0
EndEvent

% --------------------------
EventTime 5000.
Type=EndTrial
Label=end
EndEvent

OutControl

FileName=spike.all.dat
Type=Spike
population=AllPopulation
EndOutputFile

FileName=firingrate.all.dat
Type=FiringRate
population=AllPopulation
FiringRateWindow=100
PrintStep=10
EndOutputFile

FileName=firingrate.pb.dat
Type=FiringRate
population=_region_PB
FiringRateWindow=100
PrintStep=10
EndOutputFile

FileName=firingrate.ebc.dat
Type=FiringRate
population=_region_EB_C
FiringRateWindow=100
PrintStep=10
EndOutputFile
EndOutControl

