function [LevelFilt] = wellenzalThresholdFilt(Distance,Level)
%wellenzalThresholdFilt sets beat frequencies (Converted to a distance
%dimension) to Nan values when they are above the specified treshold
%
%The specified thresholds are ad hoc, and have been formed from analysing
%the Salthill dataset.
%
% Brian Scanlon, NUIG, Mar 2018

LevelFilt=Level;

%|Threshold levels:
%#1: for 0--1.5 [m], set all to Nan:
LevelFilt(Distance<1.5)=0./0;
%#2: for 1.5--3 [m], set Lvl>150 to Nan:
LevelFilt(Distance>=1.5 & Distance<3 & LevelFilt<=150)=0./0;
%#3: for 3--5.5 [m], set Lvl>80 to Nan:
LevelFilt(Distance>=3 & Distance<5.5 & LevelFilt<=80)=0./0;
%#2: for 5.5--6.5 [m], set Lvl>70 to Nan:
LevelFilt(Distance>=5.5 & Distance<6.5 & LevelFilt<=70)=0./0;
%#2: for 6.5--7.5 [m], set Lvl>60 to Nan:
LevelFilt(Distance>=6.5 & Distance<7.5 & LevelFilt<=60)=0./0;
%#2: for 7.5--8 [m], set Lvl>50 to Nan:
LevelFilt(Distance>=7.5 & Distance<8 & LevelFilt<=50)=0./0;
end

