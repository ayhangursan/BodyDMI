function [AlignedSpectra, loc]=SpectralFrequencyAlignment(Spectra,xaxis,Param,Referenceloc)
% SpectralFrequencyAlignment aligns spectra based on their frequency axis.
%
% Syntax:
%   [AlignedSpectra, loc] = SpectralFrequencyAlignment(Spectra, xaxis, Param, Referenceloc)
%
% Input Arguments:
%   - Spectra: A matrix of size [N, M] representing the input spectra, where N is the number of frequency points and M is the number of spectra.
%   - xaxis: A vector of size [N, 1] representing the frequency axis.
%   - Param: A structure containing additional parameters.
%   - Referenceloc: A vector of size [1, M] representing the reference frequency locations for alignment. If set to 0, the spectra will be aligned on a per-voxel basis.
%
% Output Arguments:
%   - AlignedSpectra: A matrix of size [N, M] representing the aligned spectra.
%   - loc: A vector of size [1, M] representing the frequency shifts applied to each spectrum for alignment.
%
% Example:
%   % Define input spectra and frequency axis
%   Spectra = rand(100, 10);
%   xaxis = linspace(0, 20, 100);
%
%   % Call the SpectralFrequencyAlignment function
%   [AlignedSpectra, loc] = SpectralFrequencyAlignment(Spectra, xaxis, Param, Referenceloc);

AlignedSpectra=zeros(size(Spectra));
% Define a range for water signal -> [-1.5 1.5] ppm
lowerend=0.2;%4.2
upperend=15.4;%15.4
water_range=find((xaxis > lowerend) & (xaxis < upperend));
disp(strcat('Water range for aligning spectra is ',num2str(lowerend),' to ',num2str(upperend),' ppm.'))
if Referenceloc==0
    disp('No reference given from perivous scan. Aligning the spectra per voxel base.' )
    
    [~, loc]=max(real(Spectra(water_range,:)),[],1);
    loc=loc+water_range(1)-1;
    Frequencyshift=xaxis(loc);
    
    for m=1:prod(Param.CSIdims)
        AlignedSpectra(:,m) = (interp1(xaxis, Spectra(:,m), xaxis-4.7+Frequencyshift(m))).';
    end
    AlignedSpectra(isnan(AlignedSpectra))=0;
    AlignedSpectra=reshape(AlignedSpectra,size(Spectra));
    loc=Frequencyshift;
else
    disp('Aligning the spectra with previous scan.' )
    
    Frequencyshift=Referenceloc;
    for m=1:prod(Param.CSIdims)
        AlignedSpectra(:,m) = (interp1(xaxis, Spectra(:,m), xaxis-4.7+Frequencyshift(m))).';
    end
    AlignedSpectra(isnan(AlignedSpectra))=0;
    AlignedSpectra=reshape(AlignedSpectra,size(Spectra));
    loc=Referenceloc;
end

end
