DatasetName=input("What is the name of the dataset you would like to Quantify(AMARES)? (ie: V1111) ","s");
if isempty(DatasetName)
    DatasetName = 'TESTdataset'; disp(['Running: ',DatasetName])
end
AMARESoptions=setAMARESoptions;
PriorKnowledge='PK_2H_7T_LiverDMI_lorentzian'; % Change to manual selection with GUI later!

eval(['fns = fieldnames(',DatasetName,');'])
if isequal(string(cellstr(fns{1})),'Noise')
    fns(1)=[];
end

filename = [DatasetName,'.xlsx'];
filenameNormalized = [DatasetName,'_Normalized.xlsx'];
ALF = 'A':'Z';
for k=1:size(fns,1)
    if isequal(strfind(fns{k},'Dyn'),1)
        eval([DatasetName,'.',fns{k},'=AMARES_masked(',DatasetName,'.',fns{k},',',DatasetName,'.DCM_Mask.DrawnROImask,AMARESoptions);'])
        MetabAmpMap(:,:,:,:,k)=eval(strcat(DatasetName,'.',string(fns(k)),'.AMARESmaps.Amplitude;'));
        for metabnum=1:size(MetabAmpMap,4)
            if k==1
                eval(['[APind, RLind, FHind]=ind2sub(',DatasetName,'.',fns{k},'.Param.CSIdims,find(',DatasetName,'.DCM_Mask.DrawnROImask>0));']);
                eval(['writecell(cellstr(strcat(num2str(APind),'','',num2str(RLind),'','',num2str(FHind))).'',''',filename,''',''Sheet'',[''Metabolite-'',num2str(metabnum)],''Range'',''A',num2str(1),''')'])%
                eval(['writecell(cellstr(strcat(num2str(APind),'','',num2str(RLind),'','',num2str(FHind))).'',''',filenameNormalized,''',''Sheet'',[''Metabolite-'',num2str(metabnum)],''Range'',''A',num2str(1),''')'])%
            end
            MetabTemp=MetabAmpMap(:,:,:,metabnum,k);
            eval(['writematrix(MetabTemp(find(',DatasetName,'.DCM_Mask.DrawnROImask>0)).'',''',filename,''',''Sheet'',[''Metabolite-'',num2str(metabnum)],''Range'',''A',num2str(k+1),''')'])%
            MetabTemp=MetabAmpMap(:,:,:,metabnum,k)./MetabAmpMap(:,:,:,end,1);
            eval(['writematrix(MetabTemp(find(',DatasetName,'.DCM_Mask.DrawnROImask>0)).'',''',filenameNormalized,''',''Sheet'',[''Metabolite-'',num2str(metabnum)],''Range'',''A',num2str(k+1),''')'])%
        end
    end
end
clear ALF MetabAmpMap MetabTemp APind RLind FHind filenameNormalized filename;
%% Save AMARES results in an Excel file

function opt=setAMARESoptions()
% Set default acquisition and experiment options
opt.PriorKnowledge='PK_2H_7T_LiverDMI_lorentzian'; % PriorKnowledge file
opt.Raw_DN=0; % Zero raw data, 1 Denoised data
prompt = {'Prior knowledge file:','Raw(0) or Denoised(1) data'};
dlgtitle = 'BodyDMI - AMARES fitting options';
fieldsize = [1 90; 1 90];
definput = {'PK_2H_7T_LiverDMI_lorentzian','0'};
answers= inputdlg(prompt,dlgtitle,fieldsize,definput);
opt.PriorKnowledge=answers{1};
opt.Raw_DN=str2double(answers{2});
end


function Dataset=AMARES_masked(Dataset,mask,AMARESoptions)
% AMARES fitting implementation for DMI in MATLAB\
% Inputs:
% inputfids: signal to be fitted in time domain
% Parameters: Acquisition parameters from DMIPipeline
% xaxis: PPM axis needed for fitting
% waterloc: Frequency location of water signal to adapt possible frequency
% offsets

% Output:
% DMIResults cell that contains amplitude, chemical shift and linewidths.

%Ayhan Gursan, UMC Utrecht 2021
% a.gursan@umcutrecht.nl
defaultfolder=cd;
%%

Parameters=Dataset.Param;
xaxis=Dataset.xaxis;
waterloc=Dataset.Referenceloc;

DMIpar.samples=Parameters.NP;
DMIpar.imagingFrequency=Parameters.Freq/(10^6);
DMIpar.timeAxis=Parameters.time;
DMIpar.dwellTime=1/(Parameters.BW);
DMIpar.ppmAxis=xaxis-4.7;
DMIpar.offset=-4.7; % This may be modified in the future to include frequency shift based on B0 imperfections
DMIpar.beginTime=Parameters.TE;disp(strcat('TE=',num2str(Parameters.TE*1000),'ms for AMARES'));

if isempty(Parameters.Index.channelIndex)
    if AMARESoptions.Raw_DN==0
        inputfids=Dataset.fftfiddata;
    elseif AMARESoptions.Raw_DN==1
        inputfids=Dataset.fftfiddataDN;
    end
else
    if AMARESoptions.Raw_DN==0
        inputfids=Dataset.RoemerComb;
    elseif AMARESoptions.Raw_DN==1
        inputfids=Dataset.RoemerCombDN;
    end
end

% Load prior knowledge file
pk=feval(AMARESoptions.PriorKnowledge);
pkglobal=pk;
CSIdims=Parameters.CSIdims;
DMI_AMARES_Results=cell(CSIdims);
for fidnum=1:numel(inputfids)/DMIpar.samples
    if isequal(mask(fidnum),1)
        for i=1:length(pk.initialValues)
            DMIpar.offset=-(4.7+(4.7-waterloc(fidnum))); % This may be modified in the future to include frequency shift based on B0 imperfections
            pk.initialValues(i).chemShift = (pkglobal.initialValues(i).chemShift + DMIpar.offset);
            pk.bounds(i).chemShift = (pkglobal.bounds(i).chemShift + DMIpar.offset);
        end
        [fitResults, fitStatus, ~, CRBResults] = AMARES.amaresFit(inputfids(:,fidnum), DMIpar, pk, 0);
        [row,col,slice]=ind2sub(CSIdims,fidnum);
        fitResults.CRLB=CRBResults;
        DMI_AMARES_Results{row,col,slice}=fitResults;
        maxfidnum=fidnum;
        cd(defaultfolder)
    end
end
disp('AMARES fitting finished.')
Dataset.DMI_AMARES_Results=DMI_AMARES_Results;
Dataset.AMARESmaps=GenerateAMARESmapsMasked(Dataset.DMI_AMARES_Results);

end