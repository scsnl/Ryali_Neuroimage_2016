clear all
close all
clc

%These scripts are Written by
% Srikanth Ryali, PhD & Vinod Menon, PhD
%Department of Psychiatry & Behavioral Sciences
%Stanford Cognitive and Systems Neuroscience Laboratory
%Stanford School of Medicine
%Stanford, USA
%Released only for Research purposes


work_directory = pwd; % Specify where the scripts are saved
addpath(genpath(work_directory))
warning('off')

%%%%%%%%%%%%%%%%%%%%% Load data Here %%%%%%%%%
%Cohort-1
%Arranged as ROI*timesamples
%ROIs - ('Motor', 'CPU', 'Insula')
%load('Data/Cohort2_Rat1.mat') %Cohort-1,Rat-1 data
load('Data/Cohort2_Rat2.mat') %Cohort-1,Rat-2 data
%%%%%%%%%%%%%%%%%%%%% Parameters Required for MDS %%%%%%%%%%
TR = 1; %For Cohort-2
tol = 10^-4;    %Tolerance for MDS convergence
maxIter = 100;   %Max allowable Iterations
%%%%%%%%%%% RUN MDS %%%%%%%%%%%%%%%%%
%To run only Motor and CPU, set timeseries = timeseries(1:2,:);
%To run CPU and Insula, set timeseries = timeseries(2:3,:);
%To Run all the three (Motor, CPU, Insula), run as it is.
main_MDS_estimation_Cohort2;
%%%%%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%%%%%%%%%%
%Zscores of the causal link
%Mapm_normal(m,n,1): Zscore for causal link from n-th ROI to m-th ROI. 
%For example:  Motor--->Thalamus, Mapm_normal(2,1,1);
Mapm_normal(:,:,1)
%P-values of the causal link
%Pvals(m,n,1): P-value for causal link from n-th ROI to m-th ROI. 
%For example:  Motor--->Thalamus, Pvals(2,1,1);
Pvals(:,:,1)
figure
plot(LL,'o-')
xlabel('Iterations')
ylabel('Log-likelihood')
title('Convergence of Log-likelihood')
figure
plot(KS.xsmooth(1,:),'linewidth',2)
hold on
plot(KS.xsmooth(2,:),'g','linewidth',2)
plot(Vm(:,1).*0.2,'r','linewidth',2)
title('Estimates of Quasi-Neuronal Signals')
legend('M1','Thalamus','Experimental Design')
