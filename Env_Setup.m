function [CPU_time_start] = Env_Setup()

    clear all; close all; clc;
    format long;
    
    % Start recording CPU time 
    CPU_time_start = clock;

    disp(' ');
    disp('********************************************************************************');
    disp('******                             GURU                                   ******');
    disp('******     General(G) Unconventional(U) Reservoir(R) Utility(U)           ******');
    disp('****         Simulator Developed by Killough Group in TAMU                ******');
    disp('****             Developer: Bicheng Yan, PhD Candidate                    ******');
    disp('****                                                                       *****');
    disp('***    Main features:                                                        ***');
    disp('***    (1) Finite Volume Method (FVM)                                        ***');
    disp('***    (2) Fully Implicit Method (FIM)                                       ***');
    disp('***    (3) Fully Compositional                                               ***');
    disp('***    (3) More physics: convection, diffusion, desorption                   ***');
    disp('***    (4) Multiple Porosity (General N Porosity)                            ***');
    disp('***    (5) Unstructured Grids / Effective Fracture Discrete Medium           ***');
    disp('********************************************************************************');    
    disp(' ');
    
end