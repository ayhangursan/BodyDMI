function patternarray=acquisitionpatterncheck(rawdata)
% acquisitionpatterncheck - Check the acquisition pattern of raw data
%
%   patternarray = acquisitionpatterncheck(rawdata) checks the acquisition
%   pattern of the given raw data and returns the pattern array. The raw
%   data should be a multi-dimensional array.
%
%   Inputs:
%       - rawdata: The raw data to check the acquisition pattern for.
%
%   Output:
%       - patternarray: The pattern array indicating the acquisition
%         pattern of the raw data.
%
%   Example:
%       rawdata = dataset.data.raw;
%       patternarray = acquisitionpatterncheck(rawdata);
%
%   Note: This function assumes that the raw data has a specific structure
%   and calculates the acquisition pattern based on the sum of absolute
%   values of the raw data along the first dimension.
dims=size(rawdata);
if numel(dims)==6
    weightedpattern=find(squeeze(sum(abs(rawdata),1)));
    B=zeros(dims(2:end));
    B(weightedpattern)=1;
    patternarray=sum(B,ndims(B));
else
    weightedpattern=find(squeeze(sum(abs(rawdata),1)));
    B=zeros(dims(2:end));
    B(weightedpattern)=1;
    patternarray=sum(B,ndims(B));    
end

% Single NSA scans result with wrong acquisition pattern matrix!
