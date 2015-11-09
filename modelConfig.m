%% Model Configuration
% Author: Xianyuan Zhan, School of Civil Engineering, Purdue University
% 
% E-mail: zhanxianyuan@purdue.edu

% This is the file that provide configuration to run the MVLPN model
% Information needs to be inputted:
% 1. Data input
% 2. Variable transformation and selection
% 3. MVPLN config
% 4. Output config
% Unless indiated, all input in this file needs to be changed in order to 
% run your MVPLN model

clear;clc;
%% Data input
% Name of the file that contains the data
% The data should be comma or tab separated txt file or Matlab .mat file, 
% with each line represent an observation and each column represent a variable 
data_file_name = 'mvarData.mat';

% Number of extra variable to be used.
% This can be used to add new variables that computed based on existing
% variables
nExtra = 1;

% Column indices of multivariate dependent variable
Y_ind = [32, 33, 34];

% Do not change following six lines
rawdata=importdata(data_file_name);
[N,M]=size(rawdata);
data=zeros(N,M+nExtra);
data(:,1:M)=rawdata;
NumSevrty = length(Y_ind);
Y=data(:,Y_ind);


%% Varialbe transformation and selection
% The existing variables or additional variables can be transformed or
% created using following method:
% Example: existing variable
% data(:,23)=data(:,23)/1000; % renormalize the scale
% data(:,12)=log(data(:,12)./(data(:,4)+data(:,5))); % linear or nonlinear transformation
% Example: additional varialbe
% data(:,M+1)=data(:,13)+data(:,14)+data(:,15); 
% Create additional variable using column index M+i, where i represents ith
% additional variable

% Column indices of the set of independent variables
selection=[0, 3, 6, 9, 10, 20, 22, 28, 30, 31]; % 0 stand for constant term
% If constant term is included, it needs to be put as the first element in the selection list


%% MVPLN config
% Maximum number of iterations
maxIter = 100;

% Number of burn-in iterations
burnin = 20;

%% Output config
% Whether save intermediate model execution data?
% These data are saved into .mat files, when loaded, can be used to
% re-compute model statistics
% 1- save; 0- not save
save_intermediate_data = 1;

% Visualization configs:
% Visualization selection choices:
% 1 for single beta variable, 2 for single sigma varialbe, 3 for all beta distribuion, 4 for sigma distribution
select=3;

% If above select choice is 1, 2 or 4,one nees to specify the variable selected
% and the corresponding depdent variable selected
% variable selected:
Var_select=9;
% Corresponding dependent variable selected
Sevrty_select=1;

% If select choice is 3, the user needs to provide the dimension for the
% subplots: dimension subpM*subpN
subpM = 3;
subpN = 3;

%% Do not change following parts
MVPLN;