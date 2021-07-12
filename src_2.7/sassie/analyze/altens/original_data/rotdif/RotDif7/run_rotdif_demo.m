%------ this is a demo script for RotDif-------
%   it loads a synthetic data set consisting of R1, R2, and NOE (here
%   called R3) and a ser of NH vectors and runs the determination of
%   the rotational diffusion tensor.
%   the relaxation data set was simulated assuming the following conditions:
%   freq=600 MHz;
%   the rotational diffusion tensor is characteruzed by the following
%   parameters: TAUc=8 ns; Dz/Dx=1.5; Dy/Dx=1.2; the Euler angles describing its
%   orientation with respect to the protein coordinate frame (pdb) are:
%   alpha=70; beta=60; gamma=170;
%   an experimental noise at the level of 2% was added to all relax. data 
%-------------------------------------------------------------------------
load rotdif_demo
RotDif(freq,R1,R2,R3,vNH);