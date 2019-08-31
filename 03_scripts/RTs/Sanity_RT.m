%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       SANITY CHECK IN RT TIMES 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Julius Welzel
% University of Oldenburg, 2018

%% Setup
clc;clear all;close all;
MAIN = ['G:\nic_LLT_all\ana03_objectbased\'];
addpath(genpath(['G:\nic_LLT_all']));
PATH_RT = [MAIN 'ana03_RTs\'];
PATH_EPOCHED = [MAIN '\ana01_epoched\'];

load ([PATH_RT 'ana02_RT_final\RT_ALL.mat']); 

%% Find all proccesed subjects

list = dir(fullfile([PATH_RT '*SUBJ*.mat']));
SUBJ = {};
for s=1:length(list)
         SUBJ = [SUBJ; list(s).name(1:6)];
end
SUBJ = unique(SUBJ);

%% Check outliers of RT with audio file

strcmp(SUBJ(1),{RT_ALL.ID})
