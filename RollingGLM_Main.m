clear; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%% This script will enable you to run analysis pipeline used in Poudel et. al. %%%
%%%   We used fMRI and behavioral data from 20 subjects exposed to rested and sleep-deprived session    %%%
%%%         If you have any queries please contact me at govinda.poudel@acu.edu.au         %%%
%%%              Below are a set of web-links needed to download external functions             %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prerequsite
% download SPM12 Toolbox

%List of subjects used in the study
Data=load('AllData.mat');

%Define the repitition time, sampling rate, and typical HRF
TR=2.5;
Fs=1/2.5;
p=[6 16 3 1 6 1 32]';
hrf=spm_hrf(2.5,p);

%Indx for ROIS on DK Atlas
task_rois=load('DK_Indx.txt');

%Spatial t-statistics map for theta activity based on Desikan Maps
thetamap=load('ThetaMap_DK.txt');
nt=length(thetamap);

% Assign data into rested and sleep-deprived session.
Subjects=Data.Data.Subjects;
Dats{1}=Data.Data.RW;
Dats{2}=Data.Data.SD;

InferredSleepEvents.Subjects=Subjects;
 InferredSleepEvents.State={'RW','SD'};

%Two sessions for each subjects
for j=1:2
    
    SubjectData=Dats{j};
    %Loop through data from all subjects
    AllState=[];
    for i=1:length(Subjects) %:5 %:length(subjects)  

            %Process RW Session Data using RollingGLM
            fmri=SubjectData.fMRI{i};
            fmri_percent=(fmri*100./mean(fmri));

            %six motion parameters
            motion=SubjectData.motion{i};
            
            %Large motion outliers generated from fsl_motion_outliers
            dvars=SubjectData.dvars{i};
            
            
            %CSF signal generated using CSF mask
            csf=SubjectData.csf{i};
            
            %Generate task regressors based on the design information
            %stored in behavioural data
            behav_data=SubjectData.behaviour{i};
            onset_volume=ceil(behav_data(:,1)/TR);
            dn=ones(length(fmri),1); % Models the mean
            dn_task=zeros(length(fmri),1);
            dn_task(onset_volume)=1;
            dn_task=conv(dn_task,hrf);
            dn_task=dn_task(1:length(fmri)); %Model task-related fMRI activity
           
            % covariates for removing from fMRI data
            covar=[dn motion dn_task dvars csf];

            % Remove the covariates from fMRI data
            resid=remove_noise(fmri_percent,covar);
            
            % Run rolling regression on the data and fit hrf at each time
            % point
            beta_values=run_rollingGLM(resid,hrf);
            
            
            % Correlate betwen beta values and EEG theta maps
            [r p]=corr(beta_values',thetamap(:));
            
            state=zeros(length(r),1);
            %Any correlation greater than 0 and p<0.05/(number of regions)
            %is a sleep-like event
            state(find(p<(0.05/nt) & r>0))=1;    
            
            AllState(:,i)=state;
            %RW_SD_States(:,i,j)=
    end
    InferredSleepEvents.Events{j}=AllState;
    
end

% A function to regress noise from fMRI signal
function bold_new=remove_noise(fMRI,regs)
    bold_new=[];
    for i=1:size(fMRI,2) 
        bold=fMRI(:,i);
        [b, int, resid]=regress(bold,regs);  
        bold_new(:,i)=resid;

    end
end

function beta_values=run_rollingGLM(fMRI,hrf)

%% Function fits the hrf model to each time point in the fMRI data in a rolling window manner
%% GLMFIT command is used to fit hrf to the same window of the data
beta_values=[];

    for i=1:size(fMRI,2)
         fMRI_data=fMRI(:,i);
         fmr=[fMRI_data(:);zeros(length(hrf),1)];
         for j=1:size(fMRI,1)

             [B1,DEV,STATS] = glmfit(hrf,fmr(j:j+length(hrf)-1));

             t=B1(2:length(B1));
             beta_values(j,i)=t(:);
         end
    end
end

