function [RoemerEqualfid, Sensitivity_maps]=RoemerEqualNoise_withmap_input(fiddata,Sensitivity_maps,NCov,channelIndex)
%
% RoemerEqualNoise_withmap_input - Roemer equal noise combination with sensitivity map as an input.
%
% Syntax: [RoemerEqualfid, Sensitivity_maps] = RoemerEqualNoise_withmap_input(fiddata, Sensitivity_maps, NCov, channelIndex)
%
% Inputs:
%   - fiddata: 5D array of complex fid data. Dimensions are [numPoints, numChannels, numLocations, numSlices, numCoils].
%   - Sensitivity_maps: 3D array of sensitivity maps. Dimensions are [numChannels, numLocations, numCoils].
%   - NCov: Noise covariance matrix.
%   - channelIndex: Index of the channel to be used for coil combination.
%
% Outputs:
%   - RoemerEqualfid: 4D array of Roemer equal noise combined fid data. Dimensions are [numPoints, AP, RL, FH].
%   - Sensitivity_maps: Updated sensitivity maps if input Sensitivity_maps is 0.
%
% Example:
%   fiddata = randn(2048, 8, 10, 10, 14);
%   Sensitivity_maps = 0;
%   NCov = eye(8);
%   channelIndex = 1;
%   [RoemerEqualfid, Sensitivity_maps] = RoemerEqualNoise_withmap_input(fiddata, Sensitivity_maps, NCov, channelIndex);
%

    if Sensitivity_maps==0
        Sensitivity_maps=(squeeze(mean((fiddata(2:5,:,:,:,:)))));disp('Sensitivity map generated') % Mean value for second to fifth points of spectrum
    else
        disp('Coil combination made with input reference sensitivity map')
    end

    dim=size(fiddata);
    grid_dims=dim(channelIndex+1:end);
    numberofloc=prod(size(Sensitivity_maps))/dim(channelIndex);
    RoemerEqualfid=zeros([dim(1) grid_dims]);
    for k=1:numberofloc
        S=Sensitivity_maps(:,k);
        U=pinv(sqrt(S'*pinv(NCov)*S))*S'*pinv(NCov);
        V=U*squeeze(fiddata(:,:,k)).';
        RoemerEqualfid(:,k)=V;
    end
end