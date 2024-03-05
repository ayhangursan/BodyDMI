function Dicom_Spectra=Show_Dicom_Spectra(Dataset)
% Show_Dicom_Spectra - Display DICOM spectra
%
%   Dicom_Spectra = Show_Dicom_Spectra(Dataset) displays DICOM spectra
%   contained in the given Dataset structure. The function plots the axial
%   image section, allows for interactive control of various parameters,
%   and shows the corresponding spectra.
%
%   Input:
%       - Dataset: A structure containing DICOM data.
%
%   Output:
%       - Dicom_Spectra: A structure containing plot properties and handles.
%
%   Example:
%       Dataset = load('dicom_data.mat');
%       Dicom_Spectra = Show_Dicom_Spectra(Dataset);
%
    fns = fieldnames(Dataset);
    DicomData=eval(strcat('Dataset.DCM_Mask;'));
    for k=1:size(fns,1)
        if ~isequal(strfind(fns{k},'Dyn'),1)
            fns{k}=[nan];
        end
    end
    fns=fns(find(cell2mat(cellfun(@(x)any(~isnan(x)),fns,'UniformOutput',false)))); % Drop noise data
    Parameters=eval(strcat('Dataset.',string(fns(1)),'.Param;'));
    inputCSI=Parameters.CSIdims;
    FIDMatrix=zeros([Parameters.NP inputCSI numel(fns)]);
    FIDMatrix_DN=zeros([Parameters.NP inputCSI numel(fns)]);
    for k=1:size(fns,1)
        if isempty(Parameters.Index.channelIndex)
            FIDMatrix(:,:,:,:,k)=eval(strcat('Dataset.',string(fns(k)),'.fftfiddata;'));
            FIDMatrix_DN(:,:,:,:,k)=eval(strcat('Dataset.',string(fns(k)),'.fftfiddataDN;'));
        else
            FIDMatrix(:,:,:,:,k)=eval(strcat('Dataset.',string(fns(k)),'.RoemerComb;'));
            FIDMatrix_DN(:,:,:,:,k)=eval(strcat('Dataset.',string(fns(k)),'.RoemerCombDN;'));
        end
    end
    FreqOffset=4.7;
    figure1=figure('WindowState','maximized');
    pause(0.4)
    figure1.Position;
    FigScaleVert=figure1.Position(4)-figure1.Position(2)+1;
    FigScaleHorz=figure1.Position(3)-figure1.Position(1)+1;

    % Axial image section
    Dicom_Spectra.PlotProp(1)=subplot(3,5,1);
    Dicom_Spectra.PlotProp(1).Position=[0.11    0.75    0.17    0.24];
    imagesc(DicomData.AxialImage(:,:,floor(inputCSI(1)/2)).^0.4)
    daspect([1 1 1])
    set(gca, 'Layer','top')
    axis on;
    [rows, columns, ~] = size(DicomData.AxialImage(:,:,1));
    hold on;
    for row = 1 : rows ./ inputCSI(1) : rows
        line([1, columns], [row, row], 'Color', 'r','Linewidth',1.5);
    end
    for col = 1 : columns ./ inputCSI(2) : columns
        line([col, col], [1, rows], 'Color', 'r','Linewidth',1.5);
    end
    colormap gray
    xtickpoints=[1:floor(columns/inputCSI(2)):columns]+floor((columns/inputCSI(2))/2);
    xticks(xtickpoints(1:2:end))
    xticklabels(string(1:2:inputCSI(2)))
    ytickpoints=[1:floor(rows/inputCSI(1)):rows]+floor((rows/inputCSI(1))/2);
    yticks(ytickpoints(1:2:end))
    yticklabels(string(1:2:inputCSI(1)))
    ax = gca;
    ax.XAxis.FontSize=16;
    ax.YAxis.FontSize=16;
    % Axial slice slider
    SliderAxial = uicontrol('style','slider','position',[FigScaleHorz*0.12 FigScaleVert*0.68 FigScaleHorz*0.14 FigScaleVert*0.02],...
        'min',1, 'max', size(DicomData.AxialImage,3),'Tag','SliderAxial', 'Value',floor(inputCSI(1)/2),'SliderStep',[1/(size(DicomData.AxialImage,3)-1) 0.1]);
    addlistener(SliderAxial, 'Value', 'PostSet', @callbackfn1);
    %% UI controls and annotations
    Apodedit = uicontrol('Style','edit','Callback',@updateeditboxvalue,'position',[FigScaleHorz*.06 FigScaleVert*0.96 60 50],'Value',0,'String',0,'FontSize',24);
    annotation(figure1,'textbox',[0.012 0.93 0.04 0.04],'String',{'Apod'},'FontSize',20,'FitBoxToText','off');
    ZFedit = uicontrol('Style','edit','Callback',@updateeditboxvalue,'position',[FigScaleHorz*.06 FigScaleVert*0.90 60 50],'Value',2,'String',2,'FontSize',24);
    annotation(figure1,'textbox',[0.012 0.87 0.04 0.04],'String',{'ZF'},'FontSize',20,'FitBoxToText','off');
    APedit = uicontrol('Style','edit','Callback',@updateeditboxvalue,'position',[FigScaleHorz*.06 FigScaleVert*0.84 60 50],'Value',floor(inputCSI(1)/2),'String',num2str(floor(inputCSI(1)/2)),'FontSize',24);
    annotation(figure1,'textbox',[0.012 0.81 0.04 0.04],'String',{'AP'},'FontSize',20,'FitBoxToText','off');
    RLedit = uicontrol('Style','edit','Callback',@updateeditboxvalue,'position',[FigScaleHorz*.06 FigScaleVert*0.83-45 60 50],'Value',floor(inputCSI(2)/2),'String',num2str(floor(inputCSI(2)/2)),'FontSize',24);
    annotation(figure1,'textbox',[0.012 0.75 0.04 0.04],'String',{'RL'},'FontSize',20,'FitBoxToText','off');
    FHedit = uicontrol('Style','edit','Callback',@updateeditboxvalue,'position',[FigScaleHorz*.06 FigScaleVert*0.81-80 60 50],'Value',floor(inputCSI(3)/2),'String',num2str(floor(inputCSI(3)/2)),'FontSize',24);
    annotation(figure1,'textbox',[0.012 0.70 0.04 0.04],'String',{'FH'},'FontSize',20,'FitBoxToText','off');

    UpdateSpectra = uicontrol('Style','pushbutton','position',[FigScaleHorz*0.03 FigScaleVert*0.65 100 50],'String','Update','FontSize',18);
    addlistener(UpdateSpectra, 'Value', 'PostSet', @callbackfn3);

    SliderPhase = uicontrol('style','slider','position',[FigScaleHorz*.005 FigScaleVert*0.55 200 50],...
        'min',-180, 'max', 180,'Tag','SliderAxial', 'Value',0,'SliderStep',[1/(2*180) 0.1]);
    addlistener(SliderPhase, 'Value', 'PostSet', @callbackfn2);
    Anotatephase=annotation(figure1,'textbox',[0.01 0.450 0.1 0.05],'String',{['Phase:',num2str(round(SliderPhase.Value))]},'FontSize',20,'FitBoxToText','off','HorizontalAlignment','center');
    Fitedit = uicontrol('Style','edit','Callback',@updateeditboxvalue,'position',[FigScaleHorz*.08 FigScaleVert*0.40 40 50],'Value',0,'String',0,'FontSize',24);
    annotation(figure1,'textbox',[0.012 0.38 0.08 0.05],'String',{'Show Fit'},'FontSize',18,'FitBoxToText','on');
    DNedit = uicontrol('Style','edit','Callback',@updateeditboxvalue,'position',[FigScaleHorz*.08 FigScaleVert*0.34 40 50],'Value',0,'String',0,'FontSize',24);
    annotation(figure1,'textbox',[0.012 0.32 0.08 0.05],'String',{'Show DN'},'FontSize',18,'FitBoxToText','on');

    % addlistener(SliderAxial, 'Value', 'PostSet', @callbackfn1);
    xaxis=eval(strcat('Dataset.',string(fns(1)),'.xaxis;')); % Raw
    show_spectra(0)

    function callbackfn1(~, eventdata)
        Axialslice = round(get(eventdata.AffectedObject, 'Value'));
        Dicom_Spectra.PlotProp(1).Children(inputCSI(1)+inputCSI(2)+1).CData = DicomData.AxialImage(:,:,Axialslice).^0.4;
        Dicom_Spectra.PlotProp(1).XLabel.String=['Axial T1w-Slice: ',num2str(Axialslice)];
    end
    function callbackfn2(~, eventdata)
        phase0 = round(get(eventdata.AffectedObject, 'Value'));
        show_spectra(phase0)
        Anotatephase.String={['Phase:',num2str(round(SliderPhase.Value))]};
    end
    function updateeditboxvalue(hObj,~)
        hObj.Value=str2double(hObj.String);
    end

    function callbackfn3(~, eventdata)
        show_spectra(SliderPhase.Value)
    end
    function show_spectra(phase0)
        ZF_factor=ZFedit.Value;
        Apodfunc=exp(-pi*Apodedit.Value.*Parameters.time).'; %Lorentzian apodization. This value could be decreased as linewidths are a bit broad now.
        if isequal(DNedit.Value,1)
            SpectraMatrix=fftshift(fft(padarray(FIDMatrix_DN.*Apodfunc,[Parameters.NP*(ZF_factor-1) 0],0,'post'),[],1),1);
        else
            SpectraMatrix=fftshift(fft(padarray(FIDMatrix.*Apodfunc,[Parameters.NP*(ZF_factor-1) 0],0,'post'),[],1),1);
        end
        if ZF_factor>1
            ppmaxis=linspace(-Parameters.ppmwindow/2,Parameters.ppmwindow/2,Parameters.NP*ZF_factor)+FreqOffset;
            FirstOrdPhase=exp(-1i * (2* pi * (ppmaxis-FreqOffset).'*(Parameters.Freq/(10^6)) * Parameters.TE));
        else
            ppmaxis=xaxis;
            FirstOrdPhase=Parameters.FirstOrdPhaseFunct;
        end

        for k=1:size(SpectraMatrix,5)
            if isequal(size(SpectraMatrix,5),1)
                subplot(3,5,[3:5 8:10 12:15]);
            else
                subplot(3,5,k+1);
            end
            cla
            Freqshift=0;
            if isequal(Fitedit.Value,1)
                SliderPhase.Value=0;%FirstOrdPhase=Parameters.FirstOrdPhaseFunct;
                tempfit=eval(strcat('Dataset.',string(fns(k)),'.DMI_AMARES_Results{APedit.Value,RLedit.Value,FHedit.Value};'));
                Freqshift=eval(strcat('Dataset.',string(fns(1)),'.Referenceloc(sub2ind(inputCSI,APedit.Value,RLedit.Value,FHedit.Value))-4.7;'));
                Resultfid=zeros(Parameters.NP,numel(tempfit.amplitude));
                time=Parameters.time+Parameters.TE;
                for metabnum=1:numel(tempfit.amplitude)
                    lineshapes(:,metabnum)=exp( -abs(time)*tempfit.linewidth(metabnum).' * pi);
                    Resultfid(:,metabnum)=(tempfit.amplitude(metabnum).*lineshapes(:,metabnum).*exp(2*pi*1i*(tempfit.chemShift(metabnum))*(Parameters.Freq/(10^6))*time).');
                end
                phase0=-tempfit.phase(1);
                FitMatrix=fftshift(fft(padarray(sum(Resultfid,2).*Apodfunc,[Parameters.NP*(ZF_factor-1) 0],0,'post'),[],1),1).*FirstOrdPhase;
                hold on
                plot(ppmaxis-Freqshift,real(FitMatrix),'b','LineWidth',2)
                plot(ppmaxis-Freqshift,real((SpectraMatrix(:,APedit.Value,RLedit.Value,FHedit.Value,k).*FirstOrdPhase.*exp(1i*pi*(phase0/180)))-FitMatrix),'r','LineWidth',2)
                hold off
                legend('Fit','Res','Raw','FontSize',8)
            end
            hold on
            plot(ppmaxis-Freqshift,real(SpectraMatrix(:,APedit.Value,RLedit.Value,FHedit.Value,k).*FirstOrdPhase.*exp(1i*pi*(phase0/180))),'k','LineWidth',2)
            hold off
            xlim([0 10])
            ylow=min(real(SpectraMatrix(:,APedit.Value,RLedit.Value,FHedit.Value,:).*FirstOrdPhase*exp(1i*pi*(phase0/180))),[],'all');
            yhigh=max(real(SpectraMatrix(:,APedit.Value,RLedit.Value,FHedit.Value,:).*FirstOrdPhase*exp(1i*pi*(phase0/180))),[],'all');
            ylim([ylow*1.2 yhigh*1.2])
            title(fns(k))
            set(gca,'XDir','reverse')
        end

    end

end