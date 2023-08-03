%% Perceptual Characterization  
clc
clear all
addpath(genpath('../'))
run loadAllDataIntoTable.m

%% Should the plots be with Weber Fraction or Perturbation Magnitude

answer = questdlg('What would you like on the x-axis?', ...
	'Options', ...
    'Weber Fraction', 'Perturbation Magnitude', 'Neither');
% Handle response
switch answer
    case 'Weber Fraction'
        disp(['Plotting against ' answer ' (|Delta V|/v)'])
        flagWF = 1;
    case 'Perturbation Magnitude'
        disp(['Plotting against ' answer ' (Delta V)'])
        flagWF = 0;    
end

%% Figure 2 and 3
close all;
[f2,f3]=accrtPlot(superSuperT, flagWF);
ph=findobj(f2,'Type','Axes');

%% Figure 4 - EZ modeling
% close all;
[f4] = EZplots(superSuperT);