function [fh,f2]=accrtPlot(trialData, flagWF)

if nargin<2 || isempty(flagWF)
    flagWF=0;  %default
end

%% Define quantities of interest

%Adding prev perturbation to table:
trialData.prevSize=[0;trialData.pertSize(1:end-1)]; 
trialData.isFirstInBlock=[1;diff(trialData.blockNo)~=0].*sign(trialData.pertSize);
trialData.prevSize(trialData.isFirstInBlock~=0)=0; %Assigning NaN to previous perturbation for first trial in each block
% trialData.prevFinalSpeed=[0;trialData.lastSpeedDiff(1:end-1)].*[0;diff(trialData.specificID)==0].*[0;diff(trialData.blockNo)==0];

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

if flagWF == 1
    %Perturbation size as a proportion of mean speed between the legs
    %(Weber Fraction)
    trialData.pertSize=trialData.pertSize./1050;
    trialData.prevSize=trialData.prevSize./1050;
end 

%% Remove first trial in each block
trialData=trialData(trialData.isFirstInBlock==0,:); 

%% Remove null and no response trials (for accurate counting of DOF in stats)
trialData=trialData(~trialData.noResponse,:);

%% Get probe sizes
B=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp=unique(trialData.pertSize); 

%% First figure: Logistic fits, accuracy and reaction times
fh=figure('Units','pixels','InnerPosition',[100 100 3*300 1*300]);
sSize=40;
[cmap,unsignedMap]=probeColorMap(23);

%% First plot Figure 1: proportion of left choices as function of probe size/weber fraction
subplot(2,3,[1,4])
hold on
set(gca,'Colormap',cmap);
S=splitapply(@(x) sum(x==-1)/sum(~isnan(x)),trialData.initialResponse,B); %Not counting NR responses
E=splitapply(@(x) nanstd(x==-1)/sqrt(sum(~isnan(x))),trialData.initialResponse,B); %Not counting NR responses
ss=scatter(pp,S,sSize,pp,'filled','MarkerEdgeColor','w');
grid on;
ylabel('proportion of left choices')
if flagWF==1
    axis([-0.360 0.360 0 1]) 
else 
    axis([-360 360 0 1]) 
end
X=trialData;
%Add fits:
hold on
errorbar(pp,S,E,'k','LineStyle','none')

X.pertSign=sign(X.pertSize);

frml='leftResponse~pertSize+isFirstInBlock+pertSize:blockNo+pertSize:pertSign+prevSize';
mm0=fitglm(X,frml,'Distribution','binomial','Link','logit')
%Automated step-down to drop non-sig terms. By default uses a deviance criterion equivalent to LRT test under Wilk's approximation
mm0=mm0.step('Upper',frml,'Criterion','Deviance','PEnter',0,'PRemove',0.05,'Nsteps',Inf)
mm0.plotPartialDependence('pertSize')
set(gca,'Colormap',unsignedMap)

%PEnter does not matter, because this is step-down strictly
%This drops the three-way pertSize:blockNo:ID interaction
%individual simple model:
Nsubs=unique(trialData.ID);
if mm0.Formula.HasIntercept
    hi='';
else
    hi='-1';
end
for i=1:length(Nsubs)
   mm{i}=fitglm(X(X.ID==Nsubs(i),:),'leftResponse~pertSize','Distribution','binomial');
   mm{i}=fitglm(X(X.ID==Nsubs(i),:),['leftResponse~' mm0.Formula.LinearPredictor hi],'Distribution','binomial');
   mm{i}.plotPartialDependence('pertSize')
end

ll=findobj(gca,'Type','Line');
set(ll(1:end-1),'Color',.7*ones(1,3));
set(ll(end),'Color','k','LineWidth',2);
uistack(ll(1:end-1),'bottom')
uistack(ss,'top')
ylabel('proportion of "left" choices') 
xlabel('R slower       same        L slower')
if flagWF==1 
    set(gca,'XLim',[-0.360 0.360])
else 
    set(gca,'XLim',[-360 360]) 
end

%% Second plot: Accuracy as a function of probe size/weber fraction
subplot(2,3,2)
set(gca,'Colormap',unsignedMap);
hold on;
grid on;
trialData.correctResponses=double(trialData.correctResponses);
trialData.correctResponses(isnan(trialData.initialResponse))=nan;
B2=findgroups(abs(trialData.pertSize));
S2=splitapply(@(x) nansum(x)/sum(~isnan(x)),trialData.correctResponses,B2); %Not counting NR responses
S2(S2==0)=NaN;
ap=sort(unique(abs(pp)));
ss=scatter(ap,S2,sSize,ap,'filled','MarkerEdgeColor','w');
E2=splitapply(@(x) nanstd(x==1)/sqrt(sum(~isnan(x))),trialData.correctResponses,B2); %Not counting NR responses
grid on
errorbar(ap,S2,E2,'k','LineStyle','none')
ylabel('accuracy') 
if flagWF==1
    axis([0 0.36 .5 1])
else
    axis([0 360 .5 1])
end

%Run binomial tests and BH on accuracy results:
clear pval h
for i=2:length(S2)
    Ntrials(i)=sum(~isnan(trialData.initialResponse(abs(trialData.pertSize)==ap(i))));
    correctTrials=Ntrials(i)*S2(i);
    pval(i-1)=binocdf(correctTrials-1,Ntrials(i),.5,'upper'); 
end
disp('Significance testing on accuracy:')
for i=1:length(pval)
    disp(['\Delta v=' num2str(ap(i+1),3) ', p=' num2str(pval(i),4) ', hits=' num2str(Ntrials(i+1)*S2(i+1),2) '/' num2str(Ntrials(i+1),2)])
end
h=BenjaminiHochberg(pval, .05); %Two-stage BKY procedure
if all(h)
    [mp,mpi]=max(pval);
    disp(['All test were significant with BKY. Largest p=' num2str(mp) ' for probe size=' num2str(ap(mpi+1)) 'mm/s'])
else
    disp('Some non-sig tests!')
    disp(num2str(ap(logical([0 h]))))
end

if flagWF==1
    xx=[0:0.001:0.360];
else
    xx=[0:360]; 
end

%Add group fit:
b=mm0.Coefficients.Estimate;
if ~mm0.Formula.HasIntercept
    y=1-1./(1+exp(b(1)*xx));
else 
   y=.5*(1-1./(1+exp(b(1)+b(2)*xx))+1./(1+exp(b(1)+b(2)*-xx)));
end
thm=find(y>.75,1,'first');
%Add individual fits:
for i=1:length(Nsubs)
    XX=X(X.ID==Nsubs(i),:);
    b=mm{i}.Coefficients.Estimate;
%     mm{i}.plotPartialDependence('pertSize')
    if ~mm{i}.Formula.HasIntercept
        y=1-1./(1+exp(b(1)*xx));
    else 
        y=.5*(1-1./(1+exp(b(1)+b(2)*xx))+1./(1+exp(b(1)+b(2)*-xx)));
    end
    
    if isempty(find(y>.75,1,'first'));
        th(i) = nan;
    else
     th(i)=find(y>.75,1,'first');
    end
end

ll=findobj(gca,'Type','Line','LineWidth',1');
uistack(ll,'bottom')
uistack(ss,'top')

disp('-------------Soft threshold stats:------------') %% I added all the xx from next line Marcela
disp(['Group=' num2str(xx(thm)) ', mean=' num2str(mean(xx(th))) ', std=' num2str(std(xx(th))) ', range=[' num2str(min(xx(th))) ',' num2str(max(xx(th))) ']']);


%% Third and forth plot: information on reaction times vs probe size/weber fraction and vs accuracy
rtPlots(trialData,[],flagWF);

%% Figure 3

k=1.0986;
X=trialData; %Excludes no response trials already
X.acc=X.correctResponses+.5*trialData.nullTrials;
biases=cellfun(@(x) x.Coefficients.Estimate(1),mm);
slopes=cellfun(@(x) x.Coefficients.Estimate(2),mm);
CI=cell2mat(cellfun(@(x) x.coefCI,mm,'UniformOutput',false));
biasCI=reshape(CI(1,:),2,9);
slopeCI=reshape(CI(2,:),2,9);
f2=figure('Units','pixels','InnerPosition',[100 100 3*300 1*300]);
subplot(1,3,2) %Bias
hold on

biases=-biases./slopes; %PSE
biasCI=-biasCI./slopes; 
gb=-mm0.Coefficients.Estimate(1)/mm0.Coefficients.Estimate(2);
gCI=-mm0.coefCI./mm0.Coefficients.Estimate(2);
if flagWF==1
    ylabel('PSE [-\beta_0/\beta_1]')
else
    ylabel('PSE [-\beta_0/\beta_1] (mm/s)')
end

%plot:
bb=bar(biases,'FaceColor',.6*ones(1,3),'EdgeColor','none');
bar(11,gb,'FaceColor',.2*ones(1,3),'EdgeColor','none');
errorbar(11,gb,gb-gCI(1,1),gCI(1,2)-gb,'k','LineStyle','none','LineWidth',1);
ee=errorbar(1:9,biases,biases-biasCI(2,:),biasCI(1,:)-biases,'k','LineStyle','none','LineWidth',1);

title('bias')
xlabel('subject')
set(gca,'XTick',[1:9,11],'XTickLabel',{'1','2','3','4','5','6','7','8','9','Group'})

subplot(1,3,3) 
hold on
%Alt: plot 1.1/beta_1
slopes=k./slopes;
slopeCI=k./slopeCI;
bb=bar(slopes,'FaceColor',.6*ones(1,3),'EdgeColor','none');
ee=errorbar(1:9,slopes,slopes-slopeCI(2,:),slopeCI(1,:)-slopes,'k','LineStyle','none','LineWidth',1);
gb=k/mm0.Coefficients.Estimate(2);
bar(11,gb,'FaceColor',.2*ones(1,3),'EdgeColor','none');
gCI=k./mm0.coefCI;
errorbar(11,gb,gb-gCI(2,2),gCI(2,1)-gb,'k','LineStyle','none','LineWidth',1);
title('probe size effect')
if flagWF==1
    ylabel('JND [1.1\beta_1^{-1}]')
else 
    ylabel('JND [1.1\beta_1^{-1}] (mm/s)')
end
set(gca,'XTick',[1:9,11],'XTickLabel',{'1','2','3','4','5','6','7','8','9','Group'})
xlabel('subject')

subplot(1,3,1) %Avg. accuracy vs. avg. RT
hold on
G=findgroups(X.ID);
acc=splitapply(@nanmean,X.acc,G);
eacc=splitapply(@(x) nanstd(x)/sqrt(numel(x)),X.acc,G);
RT=splitapply(@nanmean,X.reactionTime,G);
eRT=splitapply(@(x) nanstd(x)/sqrt(numel(x)),X.reactionTime,G);
errorbar(acc,RT,eRT,'k','LineStyle','none')
ee=errorbar(acc,RT,eacc,'k','Horizontal','LineStyle','none','DisplayName','ste');
ss=scatter(acc,RT,sSize,.5*ones(1,3),'filled','MarkerEdgeColor','w','DisplayName','Individual subject data');
rta=nanmean(X.reactionTime);
acca=nanmean(X.acc);
scatter(acca, rta,sSize,.2*ones(1,3),'filled','MarkerEdgeColor','w');
eRTa=nanstd(X.reactionTime)/sqrt(sum(~isnan(X.reactionTime)));
eacca=nanstd(X.acc)/sqrt(sum(~isnan(X.reactionTime)));
errorbar(acca,rta,eacca,'k','Horizontal','LineStyle','none','DisplayName','ste');
errorbar(acca,rta,eRTa,'k','LineStyle','none');
text(acc-.02,RT+.1,num2str([1:9]'),'Fontsize',7,'FontName','OpenSans');
text(acca+.025, rta-.2,'G','Fontsize',7,'FontName','OpenSans');
xlabel('accuracy');
ylabel('mean RT (s)');
title('indiv. accuracy vs. RT');

extendedPanelWidth(f2,.1);

hold off;

end



