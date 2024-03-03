function PhasedFIDQ2=Phase31PSpectraQ2(FID,Parameters)
% PhasedFIDQ2 = Phase31PSpectraQ2(FID, Parameters)
%
% This function applies phase correction to the input FID data using the provided parameters.
%
% Input:
%   - FID: The input FID data, a matrix of size [N, M], where N is the number of data points and M is the number of voxels, which could be in 2D or 3D.
%   - Parameters: A structure containing the phase correction parameters.
%       - FirstOrdPhaseFunct: The first-order phase correction function, a complex vector of size [N, 1].
%
% Output:
%   - PhasedFIDQ2: The phase-corrected FID data, a matrix of the same size as the input FID.
%
%   % Define the input FID data and phase correction parameters
%   FID = randn(1024, 10);
%   Parameters.FirstOrdPhaseFunct = exp(1i * pi * (0:1023)' / 1024);
%
%   % Apply phase correction
%   PhasedFIDQ2 = Phase31PSpectraQ2(FID, Parameters);

NumLines=numel(FID)./size(FID,1);
SPECTRA=reshape(fftshift(fft(FID,[],1),1),[size(FID,1) NumLines]);
FirstOrdPhaseFunct=Parameters.FirstOrdPhaseFunct;
SPECTRA=SPECTRA.*FirstOrdPhaseFunct;
PhasedSPECTRA=zeros(size(SPECTRA));
for kk=1:NumLines
    spectrum=SPECTRA(:,kk);
    phasevector=angle(spectrum);
    [~, ind]=max(abs(spectrum));
    PhasedSPECTRA(:,kk)=complex(abs(spectrum).*cos(phasevector-phasevector(ind)) ,abs(spectrum).*sin(phasevector-phasevector(ind)) );
end
PhasedSPECTRA=PhasedSPECTRA./FirstOrdPhaseFunct;
PhasedSPECTRA=reshape(PhasedSPECTRA,size(FID));
PhasedFIDQ2=ifft(ifftshift(PhasedSPECTRA,1),[],1);
