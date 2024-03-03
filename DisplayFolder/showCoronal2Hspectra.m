function showCoronal2Hspectra(inputFID,AP,ppm_axis)
% showCoronal2Hspectra displays coronal 2D spectra.
%
% Inputs:
%   - inputFID: 4D array of input FID data
%   - AP: Index of the AP (anterior-posterior) dimension
%   - ppm_axis: Vector of chemical shift values
%
% Output:
%   - None
%
% Example:
%   showCoronal2Hspectra(inputFID, 2, ppm_axis)
%
% Note:
%   - The function pads the inputFID array if its size does not match the
%     number of elements in ppm_axis.
%   - The function plots the real part of the 2D spectra.
%   - The function creates a figure with subplots for each spectrum.
%   - The x-axis represents the chemical shift values (ppm_axis).
%   - The y-axis represents the real part of the spectra.
%   - The color of the spectra is green.
%   - The background color of the subplots is black.
%   - The x-axis is reversed.
%   - The tick labels and tick marks are hidden.
%   - The x-axis limits are set to [0 10].
%   - The y-axis limits are determined based on the maximum and minimum
%     values of the spectra.

if ~isequal(size(inputFID,1),numel(ppm_axis))
    inputFID=padarray(inputFID,[(numel(ppm_axis)-size(inputFID,1)) 0 0 0],0,'post');
end

SPECTRA=fftshift(fft(permute(flip(squeeze(inputFID(:,AP,:,:)),3),[1 3 2]),[],1),1);
% Data is flipped in FH so that it will align with coronal images
% Philips CSI data slices goes from feet to head 
specmax=max(real(SPECTRA),[],'all');
% zlimitwindow=find(ppm_axis>0 & ppm_axis<10);
zlimitwindow=find(ppm_axis>-20 & ppm_axis<20);% For 31 P

specmin=1.1*min(real(SPECTRA(zlimitwindow,:,:)),[],'all');
nkx=size(SPECTRA,2);
nky=size(SPECTRA,3);

figure('color','w','WindowState','maximized')
for idx_iy = 1:nky
    for idx_ix = 1:nkx
        im_spectra = double(real(SPECTRA(:,idx_ix, idx_iy)));
        
        hAxe = axes(...
            'Parent'    , gcf                        , ...
            'Box'       , 'on'                       , ...
            'LineWidth' , 4                          , ...
            'Position'  , [(idx_iy-1)*1/nky, 1-(idx_ix*1/nkx), 1/nky, 1/nkx], ...
            'XDir'      , 'reverse'                  , ...
            'XTick'     , []                         , ...%'XTick'     , [0 2 4 6 8 10]
            'YTick'     , []                         , ...
            'XTickLabel', {''; ''}                   , ...
            'YTickLabel', {''; ''}                   , ...
            'TickDir'   , 'out'                      , ...
            'Ylim'      , [specmin specmax*1.15]     , ...
            'XLim'      , [0 10]                     ,...
            'XColor'    , [1 0 0]                    ,...
            'YColor'    , [1 0 0]                    ,...
            'Color'     , [0 0 0]);
        hLine = line( ...
            'Parent', hAxe            , ...
            'XData' , ppm_axis        , ...
            'YData' , im_spectra      , ...
            'Linewidth', 1.5         , ...
            'Color',[0 1 0]);


    end
end
end
