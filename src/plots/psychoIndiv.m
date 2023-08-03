function [xx,y,gSlope,slopes]=psychoIndiv(trialData)

%% Define quantities of interest

%Adding prev perturbation to table:
trialData.prevSize=[0;trialData.pertSize(1:end-1)]; 
trialData.isFirstInBlock=[1;diff(trialData.blockNo)~=0].*sign(trialData.pertSize);
trialData.prevSize(trialData.isFirstInBlock~=0)=0; %Assigning NaN to previous perturbation for first trial in each block
trialData.prevFinalSpeed=[0;trialData.lastSpeedDiff(1:end-1)].*[0;diff(trialData.subID)==0].*[0;diff(trialData.blockNo)==0];

%Creating binary response variable(s):
trialData.leftResponse=trialData.initialResponse==-1;
trialData.rightResponse=trialData.initialResponse==1;
trialData.noResponse=isnan(trialData.initialResponse);
trialData.nullTrials=trialData.pertSize==0;
trialData.correctResponses=trialData.initialResponse==-sign(trialData.pertSize) & ~trialData.nullTrials;
trialData.incorrectResponses=trialData.initialResponse==sign(trialData.pertSize) & ~trialData.nullTrials;

%Creating a (modified) subject ID field: %See comment at end about ANOVA
aux=trialData.subID;
trialData.ID=categorical(aux);

%Perturbation size as a proportion of mean speed between the legs
% trialData.pertSize=trialData.pertSize./1050; %Marcela
% trialData.prevSizeW=trialData.prevSize./1050; %Marcela

%% Remove first trial in each block
trialData=trialData(trialData.isFirstInBlock==0,:); 

%% Remove null and no response trials (for accurate counting of DOF in stats)
trialData=trialData(~trialData.noResponse,:);

%% Get probe sizes
B=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp=unique(trialData.pertSize); 

%% First figure: Logistic fits, accuracy and reaction times
fh=figure('Units','pixels','InnerPosition',[100 100 1*300 1*300]);
sSize=40;
[cmap,unsignedMap]=probeColorMap(23);

%% First plot Figure 1: proportion of left choices as function of probe size/weber fraction
% subplot(2,3,[1,4]);
hold on;
set(gca,'Colormap',cmap);
S=splitapply(@(x) sum(x==-1)/sum(~isnan(x)),trialData.initialResponse,B); %Not counting NR responses
E=splitapply(@(x) nanstd(x==-1)/sqrt(sum(~isnan(x))),trialData.initialResponse,B); %Not counting NR responses
ss=scatter(pp,S,sSize,pp,'filled','MarkerEdgeColor','w');
% grid on;
ylabel('proportion of left choices') 
% axis([-0.360 0.360 0 1]) %Marcela
axis([-400 400 0 1]) %Marcela
X=trialData;
%Add fits:
hold on
errorbar(pp,S,E,'k','LineStyle','none')
X.pertSign=sign(X.pertSize);

frml='leftResponse~pertSize';
mm0=fitglm(X,frml,'Distribution','binomial','Link','logit')
%Automated step-down to drop non-sig terms. By default uses a deviance criterion equivalent to LRT test under Wilk's approximation
mm0.plotPartialDependence('pertSize');
set(gca,'Colormap',unsignedMap);

xx=[-400:400];
b=mm0.Coefficients.Estimate;
if ~mm0.Formula.HasIntercept
    y{1}=1-1./(1+exp(b(1)*xx));
else
    y{1}=(1-1./(1+exp(b(1)+b(2)*xx)));
end

% plot(xx,y,'b')
gSlope=table2array(mm0.Coefficients('pertSize', 'Estimate'));

%individual simple model:
Nsubs=unique(trialData.ID);
if mm0.Formula.HasIntercept
    hi='';
else
    hi='-1';
end

slopes=[];
for i=1:length(Nsubs)
   mm{i}=fitglm(X(X.ID==Nsubs(i),:),'leftResponse~pertSize','Distribution','binomial');
   mm{i}.plotPartialDependence('pertSize')
   slopes=[slopes, table2array(mm{i}.Coefficients('pertSize', 'Estimate'))];

   b=mm{i}.Coefficients.Estimate;
   if ~mm{i}.Formula.HasIntercept
       y{i+1}=1-1./(1+exp(b(1)*xx));
   else
       y{i+1}=(1-1./(1+exp(b(1)+b(2)*xx)));
   end

end

ll=findobj(gca,'Type','Line');
set(ll(1:end-1),'Color',.7*ones(1,3));
set(ll(end),'Color','k','LineWidth',2);
uistack(ll(1:end-1),'bottom')
uistack(ss,'top')
ylabel('proportion of "left" choices') 
xlabel('R slower       same        L slower')
title('Logistic Regression')
% set(gca,'XLim',[-0.350 0.350]) %Marcela
set(gca,'XLim',[-350 350]) %Marcela



end



