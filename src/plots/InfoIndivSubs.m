%% Information on individual subjects

clc
clear all
addpath(genpath('../'))
run loadAllDataIntoTable.m

%% 
close all;
[~,~,gSlope,slopes]=psychoIndiv(superSuperT);

%%

k=gSlope*1050;

ksSE=nanstd(slopes.*1050)./sqrt(length(slopes));


% %% Is the average of the individual psychometric functions
% ... the same as the group psychometric function? 
% ytemp=[];
% for i=2:length(y)
%    
%     ytemp=[ytemp;y{i}];
%     
% end
% 
% ytemp=nanmean(ytemp);
% 
% figure();
% plot(x,y{1},'k', 'LineWidth', 2);
% hold on;
% plot(x,ytemp,'--b', 'LineWidth', 2);
% hold off;
% legend('Pooled', 'Avg');
% title('Pooled data and avergae of individual fits');
% % nanstd(1.1./slopes)./sqrt(length(slopes))