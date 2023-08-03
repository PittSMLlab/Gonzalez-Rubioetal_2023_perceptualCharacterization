function [fh]=rtPlots(trialData,goodOnly, flagWF)
if nargin<2 || isempty(goodOnly)
    goodOnly=0;
end

if nargin<3 || isempty(flagWF)
    flagWF=0;  %default
end

%%
%Creating binary response variable(s):
trialData.leftResponse=trialData.initialResponse==-1;
trialData.rightResponse=trialData.initialResponse==1;
trialData.noResponse=isnan(trialData.initialResponse);

%
trialData.absPertSize=abs(trialData.pertSize);
trialData.pertSign=sign(trialData.pertSize);

%
nullTrials=trialData.pertSize==0;
trialData.correctResponse=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; 
%
aux=trialData.subID;
aux(aux==1)=10;
aux(aux==9)=1;% To use subject 9 as reference
aux(aux==10)=9;
trialData.ID=categorical(aux);

%Define some aux vars:
correctResponses=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
incorrectResponses=trialData.initialResponse==sign(trialData.pertSize) & ~nullTrials;

B2=findgroups(abs(trialData.pertSize)); %pertSize>0 means vR>vL
pp2=unique(abs(trialData.pertSize));
Q=4;

%% Figure Features
[cmap,unsignedMap]=probeColorMap(23);
sSize=40;

%% reaction times vs. pert size
rt=trialData.reactionTime;
fun=@nanmean;

subplot(2,3,5)
RT=splitapply(fun,rt,B2); %Sort by absolute pert size
eRT=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),rt,B2); 
ss=scatter(pp2,RT,sSize,pp2,'filled','MarkerEdgeColor','w');

set(gca,'Colormap',unsignedMap)
hold on
grid on
ylabel('mean RT (s)') 
if flagWF==1
    xlabel('|\Delta V / V|')
else
    xlabel('|\Delta V| (mm/s)')
end



%Adding vR>vL and vL<vR overlapping

% set(gca,'XLim',[-10 400])
set(gca, 'XLim', [-0.1 0.381]) %Marcela
%Add" reaction times for correct vs. incorrect responses separately:
rtC=rt(correctResponses);
BC=B2(correctResponses)-1; %Subtracting 1 to eliminate the first group, which corresponds to null trials. splitapply fails for empty groups
S2=splitapply(fun,rtC,BC); %Sort by absolute pert size
E2=splitapply(@(x) std(x)/sqrt(length(x)),rtC,BC); 
rtI=rt(incorrectResponses);
BI=B2(incorrectResponses)-1;
S2=splitapply(fun,rtI,BI); %Sort by absolute pert size
E2=splitapply(@(x) std(x)/sqrt(length(x)),rtI,BI); 

e=errorbar(pp2,RT,eRT,'k');
e.LineStyle='none';
uistack(ss,'top')



%% Marcela 

% trialData.response=1-(trialData.initialResponse+1)/2;
% 
% %This comes from the ex modeling?? I ma confused
% [pSize,driftRateA,noiseA,delayA,biasA,t73A,alphaA]=fitEZ_mine(trialData);

G=findgroups(abs(trialData.pertSize));
trialData.cr=double(trialData.initialResponse==-sign(trialData.pertSize));
trialData.cr(isnan(trialData.initialResponse))=nan;
acc=splitapply(@(x) nanmean(x),trialData.cr,G);
eacc=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),trialData.cr,G);
acc(1)=.5;

rt=splitapply(@(x) nanmean(x),trialData.reactionTime,G);
ert=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),trialData.reactionTime,G);



subplot(2,3,[3,6])
ss=scatter(acc,rt,sSize,pp2,'filled','MarkerEdgeColor','w','DisplayName','group data'); %I changed pSize for pp2
set(gca,'Colormap',unsignedMap)
hold on;
grid on;
errorbar(acc,rt,ert,'k','LineStyle','none')
errorbar(acc,rt,eacc,'Horizontal','k','LineStyle','none')
hold off;
xlabel('accuracy')
ylabel('mean RT (s)')
uistack(ss,'top')
set(gca,'XLim',[.5 1.01],'YLim',[0 7]);

hold off;

%% Group fits do not work well, but individual fits may:
%Step 1: fit a psychometric curve to each subject
%Step 2: use the param fits to infer difficulty for each belt speed difference
%Step 3: use VRT, RT from all correct trials somehow to infer the best
%value for s^2 for that subject
%Idea: express VRT, MDT as a function of difficulty and (a/s). [CAN BE
%DONE], use empirical Accuracy & difficulty (or psychometric accuracy and
%difficulty?) to infer optimal value of a/s given the subject's VRT -> Does
%this work? Too little data per subject!
end

