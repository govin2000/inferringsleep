clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%% This script will enable you to run computational model of activity spread as in Poudel et. al.2020 %%%
%%%   We used fMRI and behavioral data from 20 subjects exposed to rested and sleep-deprived session    %%%
%%%         If you have any queries please contact me at govinda.poudel@acu.edu.au         %%%
%%% The model uses healthy connectome data from https://www.nitrc.org/projects/iit/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Load average connectome data from patients and controls
mat=load('Connectome.txt');

%%Altas brain ased
atlasbrain='DK_Atlas.nii';

%T-statistics of activity during sleep events in 82 rois as per Desikan
%Atlas
tstat_S1=load('RW_tstat_DK.txt'); %Rested session
tstat_S2=load('SD_tstat_DK.txt'); %Sleep deprived session


%% Set up NDM Parameters
time=0:1:20;
beta=1;
seeds=1:82;
ns=82;

%% Get Laplacian matrix
[eig_val,eig_vect, L]=generateLaplacian(mat);

%%%
%% Develop the structure for NDM
NDM.eig_vect=eig_vect;
NDM.eig_val=eig_val;
NDM.beta=beta;
NDM.time=time;
NDM.C0=zeros(ns,1);    
nVars=length(seeds);

%%%%Run NDM for all seed ROIs
[xt_all pt_all]=seedNDM(NDM,seeds);

%Repeat correaltion to generate rt curves
[rt_S1 pt]=repeatCorrelation(xt_all,tstat_S1); %Rested
[rt_S2 pt]=repeatCorrelation(xt_all,tstat_S2); %Sleep-deprived

%% Identify the time point and the region showing maximum asssocaition 
%Rested session
[rm_max_S1 im_max_S1]=max(rt_S1); 
[rmv_max_S1 imv_max_S1]=max(rm_max_S1);
tmax_S1=im_max_S1(imv_max_S1);
s_indx_s1=setdiff([1:82],imv_max_S1);

%Sleep deprived session
[rm_max_S2 im_max_S2]=max(rt_S2);
[rmv_max_S2 imv_max_S2]=max(rm_max_S2);
tmax_S2=im_max_S2(imv_max_S2);
s_indx_s2=setdiff([1:82],imv_max_S2);


[rm_sort_S1 im_sort_S1]=sort(rm_max_S1);
[rm_sort_S2 im_sort_S2]=sort(rm_max_S2);


%%Generate plot rt curve plot for maximum predictors
figure(1)
plot(rt_S1(:,imv_max_S1),'b')
hold on
plot(rt_S2(:,imv_max_S2),'r')

%%Generate scatter plot of maximum correlation
figure (2)
subplot(2,1,1)
scatter(xt_all(s_indx_s1,5,imv_max_S1),tstat_S1(s_indx_s1),'b')
subplot(2,1,2)
scatter(xt_all(s_indx_s2,5,imv_max_S2),tstat_S2(s_indx_s2),'r')

%%The scripts below generate predicted and measured activty images
nodeval_s1=rm_max_S1;
nodeval_s2=rm_max_S2;
nodeval_s1(find(nodeval_s1>-0.5))=0;
nodeval_s2(find(nodeval_s2>-0.5))=0;

%% Write predicted values for maximum association
dlmwrite('NodeVal_S1.txt',abs(nodeval_s1'),'delimiter','\t');
dlmwrite('NodeVal_S2.txt',abs(nodeval_s2'),'delimiter','\t');

%% Map of predicted activity written based on the atlas brain provided
outbrain='Predicted_activity.nii';
ros=[1:82];
xmax_S1=xt_all(:,tmax_S2,imv_max_S2)
generateFSVolumes_x(atlasbrain,ros,rescale(xmax_S1,1,10),outbrain);

%% Map of measured values as per DK atlas
outbrain='Measured_activity_S1.nii';
ros=[1:82];
generateFSVolumes_x(atlasbrain,ros,rescale(tstat_S1*1,1,10),outbrain);

outbrain='Measured_activity_S2.nii';
ros=[1:82];
generateFSVolumes_x(atlasbrain,ros,rescale(tstat_S2*1,1,10),outbrain);


%
function [xt_all pt_all]=seedNDM(NDM,seeds)
%% This function runs nework spread model for each seed region

    xt_all=[];
    pt_all=[];

    CZ=NDM.C0; %Initial condition
    for i=1:length(seeds)
        xt=[]; % 
        pt=[];
        C0=CZ;
        C0(seeds(i))=-1; %Each seed is assigned a negative value.
        NDM.C0=C0;
        [xt,pt]=RunNDM(NDM); 
        xt_all(:,:,i)=xt;
        pt_all(:,:,i)=pt;
    end

end


function [xt,pt]=RunNDM(NDM)
% This function generates network spread output for a given seed
 xt=[];
 pt=[];
 V=NDM.eig_vect;
 eig_val=NDM.eig_val;
 C0=NDM.C0;
 t=NDM.time;
 beta=NDM.beta;   
    for i=1:length(t)

       for m=1:length(eig_val)
       lambda=eig_val(m);
       u=V(:,m);   
       xst(:,m)=((exp(-beta*lambda*t(i))*u'*C0)*u);  
       pst(:,m)=(1/beta*lambda)*(1-(exp(-beta*lambda*t(i))))*(u'*C0*u);
       end
       xt(:,i)=sum(xst');
       pt(:,i)=sum(pst');

    end
end

function [eig_val V L]=generateLaplacian(mat)

%Generate laplacian matrix as per Raj et al
Deg=diag(sum(mat,2));
L=Deg-mat;

%Generate normalisation matrix 
norm_vol=generateDegreeMatrix(diag(Deg));

% normalise the laplacian matrix
L=L./norm_vol;
[V,Di]=eig(L);
eig_val=diag(Di);

end


function [VolMat]=generateDegreeMatrix(Vol)
%Finds total volume of each parcel
%
szi=size(Vol,1);
szj=size(Vol,1);

VolMat=zeros(szi,szj);

for i=1:szi
   cn1=Vol(i);
    for j=1:szj
        cn2=Vol(j);
        VolMat(i,j)=sqrt(cn1*cn2);
    end
end
end

function generateFSVolumes_x(atlasbrain,indx,atrophy,outbrain)

    V=niftiread(atlasbrain);
    infos=niftiinfo(atlasbrain);
    
    all_indx=unique(V);
    
    Vs=V;
   
    for i=1:length(indx)
        
        Vs(find(Vs==indx(i)))=atrophy(i);
        
    end
    
    intx=setdiff(all_indx,indx);
    
    for i=1:length(intx)
    
        Vs(find(Vs==intx(i)))=0;
    
    end
    niftiwrite(Vs,outbrain,infos);

end