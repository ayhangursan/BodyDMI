function outputdata=NoiseCovarianceGeneration(rawdata)
% This function generates the noise covariance matrix from the raw data.
% The raw data is read from the specified file path and processed to create
% the noise covariance matrix.

% Input:
%   - rawdata: The file path of the raw data file.

% Output:
%   - outputdata: A struct containing the processed data and noise covariance matrix.

currentfolder=cd;
outputdata.datapath=strcat(currentfolder,'\',rawdata);
[outputdata.data, outputdata.list]=csi_loadData(outputdata.datapath);

%% Parameters
outputdata.Parameters.gyromagneticratio=17.25144*10^6;
outputdata.Parameters.Tesla=7;
outputdata.Parameters.Freq=outputdata.Parameters.gyromagneticratio*outputdata.Parameters.Tesla;
outputdata.Parameters.BW=5000;
outputdata.Parameters.ppmwindow=((outputdata.Parameters.Freq/outputdata.Parameters.BW)^-1)*10^6;
outputdata.Parameters.NP=outputdata.list.F_resolution;
outputdata.Parameters.time= linspace(0,outputdata.Parameters.NP/outputdata.Parameters.BW,outputdata.Parameters.NP);
outputdata.Parameters.zerofill= outputdata.Parameters.NP*2;
outputdata.Parameters.dyn=outputdata.list.number_of_extra_attribute_1_values;

outputdata.xaxis=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.NP);

%% Concatenate of noise scan to generate a long noise scan
NoChannels=size(outputdata.data.raw,find(contains(outputdata.data.labels,'chan')));

outputdata.concatenatednoise=zeros(numel(outputdata.data.raw)/NoChannels,NoChannels);
for m=1:NoChannels
    outputdata.concatenatednoise(:,m)=reshape(squeeze(outputdata.data.raw(:,m,:)),[size(outputdata.data.raw,1)*size(outputdata.data.raw,find(contains(outputdata.data.labels,'aver'))) 1]);
end

disp('Noise scan with drive scale=0 is used to create noise covariance matrix.')
outputdata.noisecovariance=cov(outputdata.concatenatednoise);
end
