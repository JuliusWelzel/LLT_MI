%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EEG analysis of LLT group data including RTs
% Data: nic_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, Mareike Daeglau & Catharina Zich 
%(julius.welzel@uol.de, mareike.daeglau@uol.de)
% Supervisor: Cornelia Kranczioch (University of Oldenburg)

clc; clear all; close all;

MAIN = fullfile(fileparts(pwd));
addpath(genpath(MAIN));


% add toolboxes to path 
% path_toolboxes = 'C:\Users\juliu\Documents\MATLAB\toolboxes\';
% addpath([path_toolboxes 'eeglab2019_1']);
% eeglab;close all;
% addpath([path_toolboxes 'fieldtrip-20201023']);

%Change MatLab defaults
set(0,'defaultfigurecolor',[1 1 1]);

color.c_rh = [39,93,59]/255; 
color.c_lh = [117,196,107]/255;
color.c_rf = color.c_rh * 1.5;
color.c_lf = color.c_lh * 1.3;

color.c_stroke    = [104,93,121]/255;
color.c_old       = [171,108,130]/255;
color.c_young     = [216,115,127]/255;

color.c_lat = [252,187,109]/255;
color.c_med = [71,92,122]/255;



