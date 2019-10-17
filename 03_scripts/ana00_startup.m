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

MAIN = ['F:\LLT_MI\'];
addpath(genpath(MAIN));

%Change MatLab defaults
set(0,'defaultfigurecolor',[1 1 1]);

cc = cbrewer('div','RdBu',7);

c_rh = [22,162,168]/255;
c_lh = [218,40,100]/255;



