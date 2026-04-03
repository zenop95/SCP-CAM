function [pp] = switchPapCase(papCase,pp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
switch papCase
    case 1
        pp.stationKeeping          = false;  
        pp.enableSkTarget          = false;  
        pp.enableSmdGradConstraint = false;  
        pp.maxThrust               = 1e1; 
    case 2
        pp.stationKeeping          = true;  
        pp.enableSkTarget          = true;  
        pp.NCA0                    = pp.N; 
        pp.enableSmdGradConstraint = false;  
        pp.maxThrust               = 1e1; 
    case 3
        pp.stationKeeping          = true;  
        pp.enableSkTarget          = true;  
        pp.enableSmdGradConstraint = false;  
        pp.maxThrust               = 1e1; 
    case 4
        pp.stationKeeping          = false;  
        pp.enableSkTarget          = false;  
        pp.enableSmdGradConstraint = true;  
        pp.smdSoft                 = true;
        pp.maxThrust               = 1e1; 
    case 5
        pp.stationKeeping          = true;  
        pp.enableSkTarget          = true;  
        pp.enableSmdGradConstraint = true;  
        pp.maxThrust               = 10; 
        pp.smdSoft                 = false;
        pp.targWeight              = 100;                
    case 6
        pp.stationKeeping          = true;  
        pp.enableSkTarget          = true;  
        pp.enableSmdGradConstraint = true;  
        pp.smdSoft                 = true;
        pp.maxThrust               = 1e-2; 
    otherwise
end
pp.m  = 10 + pp.enableTrAndVc*12 + ...
        pp.enableSmdGradConstraint*(4 + pp.enableTrAndVc + 6*pp.smdSoft) + ...
        pp.stationKeeping*(2*pp.enableTrAndVc + 4*pp.skSoft) + ...
        pp.enableHomotopy;
pp.sl = pp.enableSkTarget*(12 + 7*pp.enableTrAndVc);
end