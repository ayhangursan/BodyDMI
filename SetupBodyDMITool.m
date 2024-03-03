% This script sets up the BodyDMI tool by performing the following steps:
% 1. Clears all variables in the workspace.
% 2. Closes all open figures.
% 3. Clears the command window.
% 4. Sets up the necessary paths for the tool.
% 5. Saves the updated path.
%
% Usage: Run this script to initialize the BodyDMI tool.

clear variables
close all
clc
%% Paths
baseExportPath = fullfile(fileparts(mfilename('fullpath')));
addpath(fullfile(baseExportPath,'ReconFolder'))
addpath(fullfile(baseExportPath,'QuantifyFolder'))
addpath(fullfile(baseExportPath,'DisplayFolder'))
addpath(fullfile(baseExportPath,'DisplayFolder/Imagescn'))
addpath(baseExportPath)
clear baseExportPath;
savepath
