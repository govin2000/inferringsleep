clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%% This script will enable you to regenerate the figures in Poudel et. al.2020 %%%
%%%   We used fMRI and behavioral data from 20 subjects exposed to rested and sleep-deprived session    %%%
%%%         If you have any queries please contact me at govinda.poudel@acu.edu.au         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prerequsite
% download SPM12 Toolbox

%Load Behavioural and Raw Data
Data=load('AllData.mat');

%Load Sleep State Data Generated by RollingGLM_main.m
Sleep=load('InferredSleep.mat');
Sleep=Sleep.InferredSleepEvents;

Subjects=Sleep.Subjects;

TR=2.5; %Repetition time of fMRI

%Define variables
ss_r_s1_rt=[];
rt_ss_s1=[];
ss_r_s2_rt=[];
rt_ss_s2=[];
avg_dur_ss_s1=[];
avg_dur_ss_s2=[];

TotalSleep_s1=sum(Sleep.Events{1}); %Total number of sleep events for each subject in rested
TotalSleep_s2=sum(Sleep.Events{2}); %Total number of sleep events for each subject in SD

for j=1:2 % Loop twice for rested and sleep deprived analysis
    ms=[];
    behav=[];
    pu=[];
    pa_ms=[];
    mean_epochs=[];

    
    if(j==1)
     ms=Sleep.Events{1};% Rested sleep events
     behavdata=Data.Data.RW.behaviour; %RW behavioural data
     pupildata=Data.Data.RW.pupil; %RW pupil data
    else
     ms=Sleep.Events{2}; %Sleep restricted sleep events
     behavdata=Data.Data.SD.behaviour; %SD behavioural data
     pupildata=Data.Data.SD.pupil; %SD pupil data   
    end
    
    mean_RT=[];
    
    for i=1:length(Subjects) %Loop through all the subjects

                
        %Run analysis for each subjects
        ss=ms(:,i);%sleep state

        %Estimate average duration of sleep events for each subject
        avg_dur_ss(i)=length(find(pulsewidth(ss)>1));

        % Find location of vols where sleep events occured
        ss_vol=find(ss>0);

        %Analyse the behavioural data 
        behav=behavdata{i};
        onset_time=behav(:,1);
        onset_volume=ceil(onset_time/TR)
        vols=zeros(length(ss),1);
        rt=behav(:,2);
        rt_vols=vols;
        rt_vols(onset_volume)=rt;
        mean_RT(i)=nanmean(rt);  %Average reaction time in response trials
        err=behav(:,3);
        total_err(i)=length(find(err==0)); %total number of errors


        diff_ss=setdiff(onset_volume,ss_vol); %Identify the trials other than sleep trials
        diff_rt(i)=nanmean(rt_vols(diff_ss)); %average RT in the trials other than sleep trials
        int_ss=intersect(onset_volume,ss_vol); %Identfy the trials coinciding with sleep trials
        int_rt(i)=nanmean(rt_vols(int_ss));% Average RT in the trials coinciding with sleep trials



        indx1=find(rt_vols>0); %Find location of trials with a response
        inx2_pa=round(ss_vol*2.5*24.5); %Find trials with sleep
        
        pu=pupildata{i};
        if(~isempty(pu))
            pa=pu(:,1);
            pa=rescale(pa,0, 1);
            pa_ms=[pa_ms; pa(find(inx2_pa<length(pa)))];
            pa_epochs=extractepochs(pa,inx2_pa,[-5,20]);
            mean_epochs(i,:)=mean(pa_epochs');

        end

    end
    
     if(j==1) %i.e. Rested
         ms_s1=ms; % Rested sleep events
         mean_RT_s1=mean_RT;
         mean_epochs_s1=mean_epochs;
         pa_ms_s1=pa_ms;
         int_rt_s1=int_rt;
         diff_rt_s1=diff_rt;
        
     else %i.e. SD
        ms_s2=ms; % SD sleep events
        mean_RT_s2=mean_RT;
        mean_epochs_s2=mean_epochs;
        pa_ms_s2=pa_ms;
        int_rt_s2=int_rt;
        diff_rt_s2=diff_rt;
    end

end

[H1,P1,CI1,STATS1] = ttest(TotalSleep_s1,TotalSleep_s2); %Paired t-test of number of sleep events between RW and SD sessions

[H2, P2, CI2, STATS2]=ttest(int_rt_s2,diff_rt_s2); %Paired t-test between RT when sleep events coincided with trials and when they didn't
[H3, P3, CI3, STATS3]=ttest(int_rt_s1,diff_rt_s1); %Paired t-test between RT when sleep events coincided with trials and when they didn't

figure(1)
subplot(2,2,1)
boxplot([TotalSleep_s1',TotalSleep_s2'])
xticklabels({'Rested','Sleep-deprived'});
ylabel('Total sleep events');

set(gca,'FontSize',18);
set(gca,'color','w');
set(gca,'box','off');
set(gcf,'color','w');
%title(['paired t-test, t=' num2str(abs(round(STATS1.tstat,2))) ' , p=' num2str(round(P1,2))]);

subplot(2,2,2)
tx=[-200:40:1000-200];
errbar=(std(mean_epochs_s1)./sqrt(18));
meanbar=mean(mean_epochs_s1);
shadedErrorBar(tx,meanbar,errbar,'lineprops','b','transparent',1);
xlabel('Time (ms)','Fontsize',18);
ylabel('Change in pupil area (%)','Fontsize',18);
axis([-200 800 -50 20])
set(gca,'FontSize',18);
set(gcf,'color','w');

hold on
errbar=(std(mean_epochs_s2)./sqrt(18));
meanbar=mean(mean_epochs_s2);
shadedErrorBar(tx,meanbar,errbar,'lineprops','r','transparent',1);
xlabel('Time (ms)','Fontsize',18);
axis([-200 800 -50 20])
set(gca,'FontSize',18);
set(gcf,'color','w');

subplot(2,2,3)
scatter(mean_RT_s1',sum(ms_s1)',200,'k','.');
lsline;
xlabel('Mean RT (ms)','Fontsize',18);
ylabel('Total sleep events','Fontsize',18);
set(gca,'FontSize',18);
set(gcf,'color','w');

subplot(2,2,4)
scatter(mean_RT_s2',sum(ms_s2)',200,'k','.');
lsline;
xlabel('Mean RT (ms)','Fontsize',18);
ylabel('Total sleep events','Fontsize',18);
set(gca,'FontSize',18);
set(gcf,'color','w');



