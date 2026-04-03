function [] = initializePath()
% initializePath Sets up MATLAB path and default graphics/UI preferences % 
% initializePath performs environment initialization used by the project: 
% - turns off the beep % - sets numeric display format to longG 
% - clears workspace 
% - adds project folders (data, Functions, CppExec, OPM, CDM) and parent 
% toolboxes to the MATLAB search path (using genpath) 
% - adds a specific external toolbox path (example: MOSEK) 
% - sets default text interpreter to LaTeX and configures default fonts, 
% font sizes, line width and figure background color
%  % % USAGE: 
% initializePath() 
% INPUT: none  
% OUTPUT: none 
% Side effects:
%
% E-mail: Author: Zeno Pavanello, 2022 
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

beep off
format longG
clear
addpath(genpath('.\data'))
addpath(genpath('.\Functions'))
addpath(genpath('.\CppExec'))
addpath(genpath('.\OPM'))
addpath(genpath('.\CDM'))
addpath(genpath('.\data_sharing'))
addpath(genpath('C:\Program Files\Mosek\10.0\toolbox\r2017a'));
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontName','Times');
set(0,'DefaultUicontrolFontName','Times', 'DefaultUicontrolFontSize', 16);
set(0,'DefaultUitableFontName','Times', 'DefaultUitableFontSize', 16);
set(0,'DefaultTextFontName','Times', 'DefaultTextFontSize', 16);
set(0,'DefaultUipanelFontName','Times', 'DefaultUipanelFontSize', 16);

set(0, 'DefaultLineLineWidth', 1);
set(0,'defaultfigurecolor',[1 1 1])
    
end