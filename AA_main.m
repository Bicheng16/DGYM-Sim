%********************************************************************************
%******                             GURU                                   ******
%******     General(G) Unstructured(U) Reservoir(R) Utility(U)             ******
%****        Simulator Developed by Killough's Group in TAMU               ******
%****             Developer: Bicheng Yan, PhD Candidate                    ******
%****                                                                       *****
%***     Main features:                                                       ***
%***    (1) Finite Volume Method (FVM)                                        ***
%***    (2) Fully Implicit Method (FIM)                                       ***
%***    (3) Fully Compositional                                               ***
%***    (3) More physics: convection(#), diffusion, desorption(#)             ***
%***    (4) Multiple Porosity (General N Porosity)                            ***
%***    (5) Unstructured Grids / Effective Fracture Discrete Medium           ***
%********************************************************************************

input = '\input\guru_5_comp.dek';

GURU(input);