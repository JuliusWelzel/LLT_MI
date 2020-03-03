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
addpath('D:\Projekt_Welzel\toolboxes\eeglab2019_0');
eeglab;close all;
addpath('D:\Projekt_Welzel\toolboxes\fieldtrip-20191111');
addpath('D:\Projekt_Welzel\toolboxes\BCILAB-devel');


%Change MatLab defaults
set(0,'defaultfigurecolor',[1 1 1]);

c_rh = [39,93,59]/255;
c_lh = [117,196,107]/255;

c_stroke    = [104,93,121]/255;
c_old       = [171,108,130]/255;
c_young     = [216,115,127]/255;

c_lat = [252,187,109]/255;
c_med = [71,92,122]/255;



