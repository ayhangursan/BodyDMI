startFolder=cd;
DatasetName=input("Welcome to Body DMI Recon script. What is the name of the dataset? (ie: V1111) ","s");
if isempty(DatasetName)
    DatasetName = 'TESTdataset'; disp(['Running: ',DatasetName])
end
NumDynamics=input("How many dynamics(including baseline) is present in the dataset? (ie: 1,2,3...)  \n");
if isempty(NumDynamics)
    disp('Number of dynamics is not given. Exiting script')
    return
elseif NumDynamics==0
    % Fill in manual raw data selection
    % NumDynamics=input("How many dynamics(including baseline) is present in the dataset? (ie: 1,2,3...)  \n input 0 to manually fill in scan numbers (ie: 15,16,17)");
elseif NumDynamics==1
    disp(['Single time point dataset'])
    pause(0.5)
    [DatasetDynName{1,1}, DatasetDynName{1,2}]=uigetfile('*.data');

elseif NumDynamics>1
    disp('Select dataset dynamics in acquisition order:')
    pause(0.5)
    DatasetDynName=cell(NumDynamics,2);
    for k=1:NumDynamics
        [DatasetDynName{k,1}, DatasetDynName{k,2}]=uigetfile('*.data');
        disp([cellstr(DatasetDynName{k,1})])
        cd(DatasetDynName{k,2})
    end
end

NumChannel=input("How many channels used during acquisition of the dataset? (ie: 1,2,3...) \n If multichannel dataset, you will select noise scan: ");
if NumChannel>1
    options=setDMIoptions(NumChannel);
    [DatasetNoiseName, DatasetNoisepath]=uigetfile('*.data');
    cd(DatasetNoisepath)
    Dataset.Noise=NoiseCovarianceGeneration(DatasetNoiseName);
    Dataset.Noise.NoiseCov=Dataset.Noise.noisecovariance;
    Dataset.Noise.NoiseCov=Dataset.Noise.NoiseCov./max(abs(Dataset.Noise.NoiseCov),[],'all');% Normalize noise covariance matrix
    options.NoiseCov= Dataset.Noise.NoiseCov;
    clear DatasetNoiseName DatasetNoisepath Dataset.Noise
elseif NumChannel==1
    disp('Single channel dataset')
    options=setDMIoptions(1);
elseif NumChannel==0
    disp('Number of coils is set to zero. Exiting script')
    return
end

%% Recons each dynamic
cd(DatasetDynName{1,2})

for k=1:NumDynamics
    eval(['Dataset.Dyn',num2str(k),'=BodyDMIRecon(''',DatasetDynName{k,1},''',options);'])
    if k==1 && NumChannel>1
        options.Referencemap=Dataset.Dyn1.RoemerSens_map;disp('First dynamic is used as baseline/Coil reference scan.')
        figure('Name','Coil array sensitivity map','WindowState','maximized')
        imagescn(abs(permute(squeeze(options.Referencemap(:,:,:,round(size(options.Referencemap,4)/2))),[2 3 1])))
        options=setDMIoptions(NumChannel, Dataset.Dyn1.options);
        if ~isequal(options.UsedCh,Dataset.Dyn1.options.UsedCh) % Check if the channels are same, if not recon again
            eval(['Dataset.Dyn',num2str(k),'=BodyDMIRecon(''',DatasetDynName{k,1},''',options);'])
            options.Referencemap=Dataset.Dyn1.RoemerSens_map;disp('First dynamic is used as baseline/Coil reference scan.')
            figure('Name','Coil array sensitivity map(Virtual channels shown. Discarded coils are set to zero)','WindowState','maximized')
            imagescn(abs(permute(squeeze(options.Referencemap(:,:,:,round(size(options.Referencemap,4)/2))),[2 3 1])))
            disp('The channels used in reconstruction have changed.')
        elseif~isequal(options.VoxelShift,Dataset.Dyn1.options.VoxelShift)
            eval(['Dataset.Dyn',num2str(k),'=BodyDMIRecon(''',DatasetDynName{k,1},''',options);'])
            disp('Voxel shift change is applied.')
        else
            disp('No changes were made in channel selection.')
        end
    end
end
clear k;
eval(([DatasetName,' = Dataset;']))
cd(startFolder)
clear Dataset DatasetName DatasetDynName NumChannel NumDynamics startFolder;


%%

function opt=setDMIoptions(NumChannel, opt)
% Set default acquisition and experiment options
if nargin < 2
    opt.TE=1.38; % Echo time
    opt.BW=2750; % Acquisition bandwidth
    opt.NoiseCov=0; % Default Noisecovariance
    opt.ReferenceLocation=0; % Default water reference signal location
    opt.Referencemap=0; % Default Roemer reference map. Zero inut generate map from the data itself.
    opt.VoxelShift=[-1 0 0];
    opt.UsedCh=[1:NumChannel].';
end
if NumChannel == 1
    prompt = {'Echo time(ms):','Acquisition bandwidth(Hz):', 'Apply Voxel shift (AP RL FH)'};
    dlgtitle = 'BodyDMI - Acquisition parameter options';
    fieldsize = [1 90; 1 90; 1 90];
    definput = {num2str(opt.TE),num2str(opt.BW),num2str(opt.VoxelShift)};
    answers= inputdlg(prompt,dlgtitle,fieldsize,definput);
    opt.TE=str2double(answers{1});
    opt.BW=str2double(answers{2});
    opt.VoxelShift=str2num(answers{3});
elseif NumChannel>1
    prompt = {'Echo time(ms):','Acquisition bandwidth(Hz):','Channels to be use in Recon (leave double space)', 'Apply Voxel shift (AP RL FH)'};
    dlgtitle = 'BodyDMI - Acquisition and Reconstruction parameter options';
    fieldsize = [1 90; 1 90; 1 90; 1 90];
    definput = {num2str(opt.TE),num2str(opt.BW),num2str(opt.UsedCh.'),num2str(opt.VoxelShift)};
    answers= inputdlg(prompt,dlgtitle,fieldsize,definput);
    opt.TE=str2double(answers{1});
    opt.BW=str2double(answers{2});
    opt.UsedCh=str2double(extract(answers{3}, digitsPattern));
    opt.VoxelShift=str2num(answers{4});
    opt.Referencemap=0; % Default Roemer reference map. Zero inut generate map from the data itself.
else
    disp("Given number of channels are not a valid number.")
end
end