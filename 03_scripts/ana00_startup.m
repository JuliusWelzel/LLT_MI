%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Setup library for Master
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EEG analysis of LLT group data including RTs
% Data: nic_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, Mareike Daeglau & Catharina Zich 
%(julius.welzel@uol.de, mareike.daeglau@uol.de)
% Supervisor: Cornelia Kranczioch (University of Oldenburg)

clc; clear all; close all;

MAIN = [fileparts(pwd) '\'];
addpath(genpath(MAIN));


% add toolboxes to path 
addpath('C:\Users\welzel-j\Desktop\toolboxes\eeglab2019_0');
eeglab;close all;
addpath('C:\Users\welzel-j\Desktop\toolboxes\fieldtrip-20191111');



%Change MatLab defaults
set(0,'defaultfigurecolor',[1 1 1]);

cc = flip(cbrewer('div','RdBu',50));

c_rh = [22,162,168]/255;
c_lh = [218,40,100]/255;



