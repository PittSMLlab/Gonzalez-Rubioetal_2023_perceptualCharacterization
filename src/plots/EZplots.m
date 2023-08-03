function [f1] = EZplots(superSuperT)

superSuperT.isFirstInBlock=[1;diff(superSuperT.blockNo)~=0].*sign(superSuperT.pertSize);
superSuperT=superSuperT(superSuperT.isFirstInBlock==0,:); 
[cmap,unsignedMap]=probeColorMap(23);

%% EZ modeling
clear driftRate
superSuperT.response=1-(superSuperT.initialResponse+1)/2;

%Group model:
[pSize,driftRateA,noiseA,delayA,biasA,t73A,alphaA]=fitEZ_mine(superSuperT);
sA=t73A/(noiseA.^2).^(1/alphaA);
aA=accFactor(pSize,biasA,sA,alphaA,noiseA);
tA=delayA+rtFactor(pSize,biasA,sA,alphaA,noiseA);

%All in a single plot:
f1=figure();
set(gcf,'Colormap',unsignedMap)
sSize=40;

%Individual model fits:
Nsubs=9;
fitMode={'RT','acc','mixed'};
for j=[2,1,3] %Different fitting modes
    clear noise aAll tAll driftRate delay bias t73 alpha
    names={'mean RT fitted', 'accuracy fitted', 'DDM fit'};
    aAll=0;
    tAll=0;
    for i=1:Nsubs
        aux=superSuperT(superSuperT.subID==i,:); %Single subj
        [pSize,driftRate(i,:),noise(i,:),delay(i),bias(i),t73(i),alpha(i)]=fitEZ_mine(aux,fitMode{j});
        a=accFactor(pSize,bias(i),t73(i),alpha(i),noise(i,:));
        aAll=aAll+a;
        t=delay(i)+rtFactor(pSize,bias(i),t73(i),alpha(i),noise(i,:)');
        tAll=tAll+t;
        difficulty(i,:)=driftRate(i,:)./noise(i,:).^2; 
       
    end
    aAll=aAll/Nsubs;
    tAll=tAll/Nsubs;
    
    if j<3
        
        subplot(2,3,4);
        hold on;
        plot(pSize,aAll,'LineWidth',2,'DisplayName', names{j});
        
        subplot(2,3,5);
        hold on;
        plot(pSize,tAll,'LineWidth',2,'DisplayName',names{j});
        
        subplot(2,3,6);
        hold on;
        plot(aAll,tAll,'LineWidth',2,'DisplayName',names{j});
        
    else
        
        subplot(2,3,1);
        hold on;
        p1 = plot(pSize,aAll,'k','LineWidth',2,'DisplayName',names{j});
        subplot(2,3,2);
        hold on;
        plot(pSize,tAll,'k','LineWidth',2,'DisplayName',names{j});
        
        subplot(2,3,3);
        hold on;
        plot(aAll,tAll,'k','LineWidth',2,'DisplayName',names{j});
        
    end
end

%Then plot data:
subplot(2,3,1);
hold on;
G=findgroups(abs(superSuperT.pertSize));
superSuperT.cr=double(superSuperT.initialResponse==-sign(superSuperT.pertSize));
superSuperT.cr(isnan(superSuperT.initialResponse))=nan;
acc=splitapply(@(x) nanmean(x),superSuperT.cr,G);
eacc=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),superSuperT.cr,G);
acc(1)=.5;
ss=scatter(pSize,acc,sSize,pSize,'filled','MarkerEdgeColor','w'); %all data
errorbar(pSize,acc,eacc,'k','LineStyle','none');
xlabel('|DeltaV/V| (mm/s)');
ylabel('accuracy');
set(gca,'YLim',[.5 1.01]);
uistack(ss,'top');

subplot(2,3,4);
hold on;
G=findgroups(abs(superSuperT.pertSize));
superSuperT.cr=double(superSuperT.initialResponse==-sign(superSuperT.pertSize));
superSuperT.cr(isnan(superSuperT.initialResponse))=nan;
acc=splitapply(@(x) nanmean(x),superSuperT.cr,G);
eacc=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),superSuperT.cr,G);
acc(1)=.5;
ss=scatter(pSize,acc,sSize,pSize,'filled','MarkerEdgeColor','w'); %all data
errorbar(pSize,acc,eacc,'k','LineStyle','none');
xlabel('|DeltaV/V| (mm/s)');
ylabel('accuracy');
set(gca,'YLim',[.5 1.01]);
uistack(ss,'top');

subplot(2,3,2);
hold on;
rt=splitapply(@(x) nanmean(x),superSuperT.reactionTime,G);
ert=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),superSuperT.reactionTime,G);
ss=scatter(pSize,rt,sSize,pSize,'filled','MarkerEdgeColor','w'); %all data
errorbar(pSize,rt,ert,'k','LineStyle','none');
xlabel('|\Delta V| (mm/s)');
ylabel('mean RT (s)');
legend('Location','NorthEast','Box','off');
set(gca,'YLim',[0 7]);
uistack(ss,'top');

subplot(2,3,5);
hold on;
rt=splitapply(@(x) nanmean(x),superSuperT.reactionTime,G);
ert=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),superSuperT.reactionTime,G);
ss=scatter(pSize,rt,sSize,pSize,'filled','MarkerEdgeColor','w'); %all data
errorbar(pSize,rt,ert,'k','LineStyle','none');
xlabel('|\Delta V| (mm/s)');
ylabel('mean RT (s)');
legend('Location','NorthEast','Box','off');
set(gca,'YLim',[0 7]);
uistack(ss,'top');

subplot(2,3,3);
hold on;
ss=scatter(acc,rt,sSize,pSize,'filled','MarkerEdgeColor','w'); %all data
errorbar(acc,rt,ert,'k','LineStyle','none');
errorbar(acc,rt,eacc,'Horizontal','k','LineStyle','none');
xlabel('accuracy');
ylabel('mean RT (s)');
uistack(ss,'top');
set(gca,'XLim',[.5 1.01],'YLim',[0 7]);


subplot(2,3,6);
hold on;
ss=scatter(acc,rt,sSize,pSize,'filled','MarkerEdgeColor','w'); %all data
errorbar(acc,rt,ert,'k','LineStyle','none');
errorbar(acc,rt,eacc,'Horizontal','k','LineStyle','none');
xlabel('accuracy');
ylabel('mean RT (s)');
uistack(ss,'top');
set(gca,'XLim',[.5 1.01],'YLim',[0 7]);

hold off;

end

