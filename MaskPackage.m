DatasetName=input("What is the name of the dataset you would like to display? (ie: V1111) ","s");
% This script prompts the user to enter the name of a dataset and then
% loads the DICOM images and draws a mask on them. The script also handles
% the case when the user does not provide a dataset name by setting a
% default name and displaying a message. After loading and drawing the
% mask, the script waits for the user to close the figure window before
% continuing. Finally, the script assigns the mask to the dataset's
% DCM_Mask property and clears unnecessary variables.
%
% Inputs:
%   - None
%
% Outputs:
%   - None
%
% Usage:
%   1. Run the script.
%   2. Enter the name of the dataset when prompted.
%   3. Close the figure window to continue.
%
% Example:
%   What is the name of the dataset you would like to display? (ie: V1111) V1234
%
%   Running: V1234
%
%   (Figure window opens with DICOM images and mask drawn)
%
%   (User closes the figure window)
%
%   (Script continues executing)

if isempty(DatasetName)
    DatasetName = 'TESTdataset'; disp(['Running: ',DatasetName])
end
eval(([DatasetName,'.DCM_Mask=LoadDICOMandDrawMask(',DatasetName,');']))
fig = gca;
if isvalid(fig)
    waitfor(fig);
end
eval(([DatasetName,'.DCM_Mask=DCM_Mask;']))
clear DCM_Mask fig