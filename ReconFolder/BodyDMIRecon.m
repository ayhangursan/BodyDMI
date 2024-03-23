function dataset=BodyDMIRecon(rawdata,options)
% Pipeline for 31P CSI experiments with multichannel acquisition
% made in Philips scanner(7T)
% Example use V9999=P31Pipeline('raw_012','raw_013',0.75,5000)
%'raw_012' is CSI data, 'raw_013' is seperate noise scan
% 0.75 is TE in ms, 5000 is BW in Hz

% Outline
% 1-Read dataset
% 2-Define parameters
% 3-Average data
% 4-Apply additional Hamming weighting correction(already k-space weighted acqusition)
% 5-Spatial FFT
% 6-Voxel base phasing spectra and channel combination
% 7-Spectral FFT(if it is not already applied)
% 8-Denoising on one dataset(Could be modified)
% Notes
%% Ayhan Gursan 2023, agursan@umcutrect.nl

dataset.options = options;
disp(strcat('Processing CSI data:',rawdata));

% TE=1.49;% Echo time is hard coded
%Read seperate noise scan
dataset.NoiseCov=options.NoiseCov;
%Read Philips .data file
currentfolder=cd;
dataset.datapath=strcat(currentfolder,'\',rawdata);
[dataset.data, dataset.list]=csi_loadData(dataset.datapath);
%% Parameters
dataset.Param.gyromagneticratio=6.53569*10^6;
dataset.Param.Tesla=7;
dataset.Param.Freq=dataset.Param.gyromagneticratio*dataset.Param.Tesla;
if nargin < 5
    dataset.Param.BW=2750;disp(strcat('Default spectral BW:',num2str(dataset.Param.BW)));
else
    dataset.Param.BW=options.BW;disp(strcat('Spectral BW:',num2str(dataset.Param.BW)));
end
%%
dataset.Param.ppmwindow=((dataset.Param.Freq/dataset.Param.BW)^-1)*10^6;
dataset.Param.NP=dataset.list.F_resolution;
dataset.Param.time= linspace(0,dataset.Param.NP/dataset.Param.BW,dataset.Param.NP);
ZerofillFactor=8;
dataset.Param.zerofill= dataset.Param.NP*ZerofillFactor;

FreqOffset=4.7;
dataset.xaxis=linspace(-dataset.Param.ppmwindow/2,dataset.Param.ppmwindow/2,dataset.Param.NP)+FreqOffset;
dataset.xaxiszerofill=linspace(-dataset.Param.ppmwindow/2,dataset.Param.ppmwindow/2,dataset.Param.zerofill)+FreqOffset;
dataset.Param.TE=options.TE*10^-3; % Echo time as an input. Reverse it later !!!
disp(strcat('Acquisition echo time:', num2str(dataset.Param.TE*10^3),' ms'))
dataset.Param.missingpoints=ceil(dataset.Param.BW*dataset.Param.TE);
dataset.Param.Apodizationparameter=5; % Line broadening parameter(in Hertz)
dataset.Param.apodfunc=exp(-pi*dataset.Param.Apodizationparameter.*dataset.Param.time).'; %Lorentzian apodization. This value could be decreased as linewidths are a bit broad now.
% Index
dataset.Param.Index.channelIndex = find(contains(dataset.data.labels,'chan'));
dataset.Param.Index.averageIndex = find(contains(dataset.data.labels,'aver'));
dataset.Param.Index.kxIndex = find(contains(dataset.data.labels,'kx'));
dataset.Param.Index.kyIndex = find(contains(dataset.data.labels,'ky'));
dataset.Param.Index.kzIndex = find(contains(dataset.data.labels,'kz'));
dataset.Param.FirstOrdPhaseFunct=exp(-1i * (2* pi * (dataset.xaxis-FreqOffset).'*(dataset.Param.Freq/(10^6)) * dataset.Param.TE));
dataset.Param.FirstOrdPhaseFunctZF=exp(-1i * (2* pi * (dataset.xaxiszerofill-FreqOffset).'*(dataset.Param.Freq/(10^6)) * dataset.Param.TE));
%% Control if data size and parameters in list file match

if ~isequal(size(dataset.data.raw,dataset.Param.Index.kxIndex),dataset.list.X_resolution)
    disp('Number of encodings in raw data does not match with acquisition. Probably shutter was used in this acquisition. Zeropadding the data in X axis')
    Difference_X=dataset.list.X_resolution-size(dataset.data.raw,dataset.Param.Index.kxIndex);
    tempindex=zeros(numel(size(dataset.data.raw)),1);
    tempindex(dataset.Param.Index.kxIndex)=1;
    dataset.data.raw=padarray(dataset.data.raw,tempindex*ceil(Difference_X/2),0,'pre');
    dataset.data.raw=padarray(dataset.data.raw,tempindex*floor(Difference_X/2),0,'post');
end
if ~isequal(size(dataset.data.raw,dataset.Param.Index.kyIndex),dataset.list.Y_resolution)
    disp('Number of encodings in raw data does not match with acquisition. Probably shutter was used in this acquisition. Zeropadding the data in Y axis')
    Difference_Y=dataset.list.Y_resolution-size(dataset.data.raw,dataset.Param.Index.kyIndex);
    tempindex=zeros(numel(size(dataset.data.raw)),1);
    tempindex(dataset.Param.Index.kyIndex)=1;
    dataset.data.raw=padarray(dataset.data.raw,tempindex*ceil(Difference_Y/2),0,'pre');
    dataset.data.raw=padarray(dataset.data.raw,tempindex*floor(Difference_Y/2),0,'post');
end
if ~isequal(size(dataset.data.raw,dataset.Param.Index.kzIndex),dataset.list.Z_resolution)
    disp('Number of encodings in raw data does not match with acquisition. Probably shutter was used in this acquisition. Zeropadding the data in Z axis')
    Difference_Z=dataset.list.Z_resolution-size(dataset.data.raw,dataset.Param.Index.kzIndex);
    tempindex=zeros(numel(size(dataset.data.raw)),1);
    tempindex(dataset.Param.Index.kzIndex)=1;
    dataset.data.raw=padarray(dataset.data.raw,tempindex*ceil(Difference_Z/2),0,'pre');
    dataset.data.raw=padarray(dataset.data.raw,tempindex*floor(Difference_Z/2),0,'post');
end

dataset.Param.dims=size(dataset.data.raw);
dataset.Param.CSIdims=[dataset.Param.dims(dataset.Param.Index.kxIndex) dataset.Param.dims(dataset.Param.Index.kyIndex) dataset.Param.dims(dataset.Param.Index.kzIndex)];
disp(['Eventual data MRSI grid X:',num2str(dataset.Param.CSIdims(1)),' Y:',num2str(dataset.Param.CSIdims(2)),' Z:',num2str(dataset.Param.CSIdims(3))])

%% Average data
if dataset.list.number_of_signal_averages==1
    dataset.avgrawdata=squeeze((dataset.data.raw)); disp('NSA=1')
    % Spatial filter-Hamming for 3D
    dataset.acqpattern=acquisitionpatterncheck(dataset.data.raw);
    dataset.acqpattern=ones(size(dataset.acqpattern));
else
    dataset.avgrawdata=squeeze(mean(dataset.data.raw,dataset.Param.Index.averageIndex));
    % Spatial filter-Hamming for 3D
    dataset.acqpattern=acquisitionpatterncheck(dataset.data.raw);

end


if isempty(dataset.Param.Index.channelIndex)==0
    dataset.acqpattern=squeeze(dataset.acqpattern(1,:,:,:)/max(dataset.acqpattern,[],'all'));
end

[dataset.idealhammingwindow,dataset.Correctionmask,dataset.avgrawdata]=Hammingfilterpostcorrection(dataset.avgrawdata,dataset.acqpattern,dataset.Param);
% FFT for each channel
% dataset.fftfiddata=zeros(size(dataset.avgrawdata));
% dataset.spectradata=zeros(size(dataset.avgrawdata)); % Save memory
if ~isempty(dataset.Param.Index.channelIndex)
    dataset.avgrawdata=permute(dataset.avgrawdata,[1 2 4 3 5])*1e3;disp('Signal rescaled with 1e3, fix amount!')
    dataset.avgrawdata=dataset.avgrawdata(:,options.UsedCh(:),:,:,:);disp('Removed excluded channels, if any.')
    dataset.Param.dims(2)=numel(options.UsedCh);
    dataset.Param.CSIdims=[dataset.Param.dims(dataset.Param.Index.kyIndex) dataset.Param.dims(dataset.Param.Index.kxIndex) dataset.Param.dims(dataset.Param.Index.kzIndex)];
    dataset.avgrawdata=flip(dataset.avgrawdata,3);disp('Fliped in AP')
    dataset.avgrawdata=flip(dataset.avgrawdata,4);disp('Fliped in RL')
    dataset.fftfiddata=zeros(size(dataset.avgrawdata));
    for n=1:size(dataset.avgrawdata,dataset.Param.Index.channelIndex)
        dataset.fftfiddata(:,n,:,:,:)=circshift(PhaseSpectra(SpatialFFT(squeeze(dataset.avgrawdata(:,n,:,:,:))),dataset.Param),[0 options.VoxelShift]);
    end
    disp(['Applied voxel shift(AP RL FH) = ',num2str(options.VoxelShift)])
    disp(['Channels used in reconstruction: ',num2str(options.UsedCh.')])
else
    dataset.avgrawdata=permute(dataset.avgrawdata,[1 3 2 4])*1e3;disp('Signal rescaled with 1e3, fix amount!')
    dataset.Param.CSIdims=[dataset.Param.dims(dataset.Param.Index.kyIndex) dataset.Param.dims(dataset.Param.Index.kxIndex) dataset.Param.dims(dataset.Param.Index.kzIndex)];
    dataset.avgrawdata=flip(dataset.avgrawdata,2);disp('Fliped in AP')
    dataset.avgrawdata=flip(dataset.avgrawdata,3);disp('Fliped in RL')
    dataset.fftfiddata=circshift(PhaseSpectra(SpatialFFT(squeeze(dataset.avgrawdata)),dataset.Param),[0 options.VoxelShift]);
    dataset.fftfiddata=Phase31PSpectraQ2(dataset.fftfiddata,dataset.Param);
    disp(['Applied voxel shift(AP RL FH) = ',num2str(options.VoxelShift)])
end
disp('Spatial FFT applied.')

    removedfield='avgrawdata';dataset=rmfield(dataset,removedfield);disp('Removing avgrawdata for clearing memory.')
    removedfield='data'; dataset=rmfield(dataset,removedfield);disp('Removing RAW data for clearing memory.')
%%
if ~isempty(dataset.Param.Index.channelIndex) % if multichannel dataset
    [dataset.DecorrelatedSignal, ~]=DecorrDenoise(dataset.fftfiddata,dataset.NoiseCov,dataset.Param);
    % Phase each spectra
    dataset.PhasedFID=Phase31PSpectraQ2(dataset.DecorrelatedSignal,dataset.Param);

    disp('Automatic phasing applied.')
    % Channel combination
    [dataset.RoemerComb, dataset.RoemerSens_map]=RoemerEqualNoise_withmap_input(dataset.PhasedFID,options.Referencemap,diag(ones(dataset.Param.dims(dataset.Param.Index.channelIndex),1)),dataset.Param.Index.channelIndex);        disp('Roemer equeal noise channel combination applied.')
    dataset.RoemerComb=PhaseSpectra(dataset.RoemerComb,dataset.Param); % Phase combined data
%     dataset.Combspectra=fftshift(fft(PhaseSpectra(dataset.RoemerComb,dataset.Param),[],1),1).*dataset.Param.FirstOrdPhaseFunct; % Phase and apply spectral FFT Roemer
    dataset.RoemerCombDN=PCA_CSIdenoising_V2_KM(dataset.RoemerComb,5,dataset.Param); % First combined and then denoised.

%     [dataset.DN_RoemerComb, dataset.DNRoemerSens_map]=RoemerEqualNoise_withmap_input(dataset.DNPhasedFID,0,diag(ones(dataset.Param.dims(dataset.Param.Index.channelIndex),1)),dataset.Param.Index.channelIndex);        disp('Roemer equeal noise channel combination applied.')
%     dataset.DNCombspectra=fftshift(fft(Phase31PSpectraQ2(dataset.DN_RoemerComb,dataset.Param),[],1),1).*dataset.Param.FirstOrdPhaseFunct; % Phase and apply spectral FFT Roemer
    removedfield='fftfiddata'; dataset=rmfield(dataset,removedfield);disp('Removing fftfiddata(multichannel before comb) data for clearing memory.')
    removedfield='DecorrelatedSignal'; dataset=rmfield(dataset,removedfield);disp('Removing DecorrelatedSignal data for clearing memory.')
    removedfield='PhasedFID'; dataset=rmfield(dataset,removedfield);disp('Removing PhasedFID(before comb) data for clearing memory.')

    %% SNR estimation on raw data(no denoising)
    % Currently SNR is not taken into account in any process, however it
    % can still be calculated once and saved without occupiying memory at a
    % significant size.
    WaterFreq=0;
    dataset.noisewindow=find(dataset.xaxis+WaterFreq>15 & dataset.xaxis+WaterFreq<25);
    dataset.waterwindow=find(dataset.xaxis+WaterFreq>3.8 & dataset.xaxis+WaterFreq<5.7);
    temp=fftshift(fft(dataset.RoemerComb,[],1),1).*dataset.Param.FirstOrdPhaseFunct; % Apply spectral FFT to calculate SNR
    dataset.noisemap=squeeze(std(real(temp(dataset.noisewindow,:,:,:)))); % Using real part of the spectra to estimate noise in accordance with MRS consensus paper(doi:10.1002/nbm.4347)i
    dataset.SNR=squeeze(max(real(temp(dataset.waterwindow,:,:,:)),[],1))./dataset.noisemap;

    [~, dataset.Referenceloc]=SpectralFrequencyAlignment(abs(temp),dataset.xaxis,dataset.Param,options.ReferenceLocation);
else
    dataset.fftfiddataDN=PCA_CSIdenoising_V2_KM(dataset.fftfiddata,5,dataset.Param);
    WaterFreq=0;
    dataset.noisewindow=find(dataset.xaxis+WaterFreq>15 & dataset.xaxis+WaterFreq<25);
    dataset.waterwindow=find(dataset.xaxis+WaterFreq>3.8 & dataset.xaxis+WaterFreq<5.7);
    temp=fftshift(fft(dataset.fftfiddata,[],1),1).*dataset.Param.FirstOrdPhaseFunct; % Apply spectral FFT to calculate SNR
    dataset.noisemap=squeeze(std(real(temp(dataset.noisewindow,:,:,:)))); % Using real part of the spectra to estimate noise in accordance with MRS consensus paper(doi:10.1002/nbm.4347)i
    dataset.SNR=squeeze(max(real(temp(dataset.waterwindow,:,:,:)),[],1))./dataset.noisemap;
    [~, dataset.Referenceloc]=SpectralFrequencyAlignment(real(temp),dataset.xaxis,dataset.Param,options.ReferenceLocation);
end
end

function [DecorrelatedSignal, DenoisedDecorrelatedSignal]=DecorrDenoise(Signal,NoiseCov,Param)
% Compact function to decorrelate signal and denoise afterwards.
% Decorrelation is made with Cholesky decomposition of noise covariance
% matrix.
% Denoising is done with PCA-Denoising.

% Ayhan Gursan, 2021, UMC Utrecht.
% a.gursan@umcutrecht.nl


% 1- Use Cholesky decomposition for decorrelation
% CG-SENSE, Pruessman et al MRM 2001.
% The basic idea is to create a set of virtual channels by linear
% combination of the original ones, such that the virtual channels exhibit
% unit noise levels and no mutual noise correlation.

% We will use Cholesky decompostion to generate these virtual coils that
% decorrelated and have unit noise.
L=chol(NoiseCov);
% Matlab gives the hermetian of the decomposition as output so in the following steps L' will be used for decorrelation of signal.
% For validation check NoiseCov=L'*L;
dims=size(Signal);
channelIndex=Param.Index.channelIndex;

DecorrelatedSignal=reshape(Signal,[dims(1) dims(channelIndex) numel(Signal)/(dims(1)*dims(channelIndex))]);
for kk=1:numel(Signal)/(Param.NP*Param.dims(2))
    DecorrelatedSignal(:,:,kk)=squeeze(pinv(L')*DecorrelatedSignal(:,:,kk).').';
end
DecorrelatedSignal=reshape(DecorrelatedSignal,dims);

% 2- Denoise each channel with PCA-Denoising
% PCA-Denoising, Froeling et al MRM 2021.
ReconFast=1;
if ReconFast==1
    DenoisedDecorrelatedSignal=ones(size(DecorrelatedSignal));disp('No per channel denoising applied.')
    DenoisedDecorrelatedSignal=DenoisedDecorrelatedSignal.*rand(size(DenoisedDecorrelatedSignal));
else
    DenoisedDecorrelatedSignal=zeros(size(DecorrelatedSignal));
    patchsize=5;
    for ch=1:dims(2)
        DenoisedDecorrelatedSignal(:,ch,:,:,:) = PCA_CSIdenoising_V2_KM(squeeze(DecorrelatedSignal(:,ch,:,:,:)),patchsize,Param);
    end

    % DenoisedDecorrelatedSignal could be used in channel combination with
    % identity matrix instead of noise covariance matrix.

    %  Such as
    % [Decorr.Roemer, Decorr.RoemerSensitivityMaps] = RoemerEqualNoise_withmap_input(DenoisedDecorrelatedSignal,0,diag(ones(16,1)),Param.Index.channelIndex);% Roemer combination with decorrelated DN signal
    % or
    % [Decorr.WSVD, Decorr.WSVDQuality, Decorr.WSVDWeights]=WSVDcomb(DenoisedDecorrelatedSignal,diag(ones(16,1)),Param.Index.channelIndex);        % WSVD combination with decorrelated DN signal
end
end

function PhasedFID=PhaseSpectra(FID,Parameters)
% Phase signal in frequency domain by matching the real with magnitude at
% maximum intensity point.

% Input is time domain signal!

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
PhasedFID=ifft(ifftshift(PhasedSPECTRA,1),[],1);
end
