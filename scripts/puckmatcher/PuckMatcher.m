function PuckMatcher(scriptfolder, manifest, mypuck)
%%%%This pipeline function assumes that we are using the large image
%%%%feature and the XY feature in Nikon images, and that files ending in
%%%%xy1 are for puck 1. We do not do stitching.

%This function is for beads with 7 J bases in two groups.

%%%%SETUP:

%0) Make sure ImageJ is closed.

%1) When you are done with the code, you should make the raw data folder in
%Pucks and the InputFolder in find_roi online only via smart sync to save
%hard drive space. You should also delete the pucktmp directory.

%2) A number of paths are hardcoded in the code currently, e.g.
%"C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code\find_roi\helpers\vlfeat-0.9.20\toolbox\"
%in find_roi_stack_fun

%3) Change ImageSize to reflect something slightly smaller than the size of
%the final stitched images.

%4) If you have a digital gene expression matrix, you need to take the entire CSV and put it in the Illumina folder (C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\Pucks\Illumina\Puck_180106_3) and call it DGE.csv. I then copy out the top row with notepad++, delete the first entry (which is there as a placeholder for the row labels) and put it into a different file called IlluminaBarcodes.txt

%5) Each ligation sequence and primer set will require a bunch of modifications
%to MapLocationsFun. To accommodate for this, I make a different version of MapLocationsFun for
%each ligation sequence. You have to change which version is called below.

%6) NOTE: SlideseqSave.ijm is stored in C:\Fiji.app\macros, and the text
%is:
%params=getArgument()
%print("Opening file:")
%print(params)
%open(params);
%saveAs("Tiff","C:\\PuckOutputTmp.tif")
%print("Done.")
%exit();


%For future releases, consider setting the priority of the parpool
%processes to High:

%matlabpool open
%cmd_str = 'wmic process where name="MATLAB.exe" CALL setpriority 64';
%[~,~] = system(cmd_str);


%% Initialize
%Set up the paths to all of the various folders with functions in them that
%we will call:
%clear all
%close all

%addpath([scriptfolder, '/puckmatcher']);
%PythonPath='C:\Users\sgr\AppData\Local\Programs\Python\Python37\python.exe';
%BeadseqCodePath='D:\Jilong\SlideseqCode\BeadSeq Code';
%PipelineFunctionPath='D:\Jilong\SlideseqCode\PipelineFunctions';
%addpath(BeadseqCodePath,[BeadseqCodePath,'\find_roi'],PipelineFunctionPath);
%addpath([BeadseqCodePath,'\find_roi\helpers']);
%We assume that the nd2s have been exported to tiffs in the format:
%DescriptiveNametLYX, where DescriptiveName refers to one run of the microscope,  Y is a letter and X is a number, and L is the
%name of the ligation within that DescriptiveName file series.
%We convert the files to a final output format:
%Puck85 Ligation X Position AB
%BeadType="180402"; for 14bp barcodes from 180402 beads
BeadType="180728";

RenameFiles=1;

RunSignificanceAnalysis=1;  %not important
NumClusters=[11,17,14,14]; %Cerebellum is 1, hippocampus is 2, frontalcortex is 3, posteriorcortex is 4 % not important 

%Which pucks do we try to run analogizer on? 0 if don't run, otherwise, the
%number here refers to the element of DropseqDGEPaths and
%DropseqClusterPaths to use as a reference.
RunAnalogizer=[0,2,2,0,2,0,2,2,2,2,2,2,2,0,0,0,1,1,1,0,2,2,2,2]; %not i,portant
AnalogizerBeadCutoff=5; %Note that this threshold is only applied for the variable genes used in the NMFreg, so it can be low
AnalogizerType="NMFReg";
%Cerebellum is 1. Hippocampus is 2.
%DropseqDGEPaths={'\\helium\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Cerebellum_ALT','\\iodine-cifs\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Hippocampus'};
%DropseqClusterPaths={'\\helium\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Cerebellum_ALT\assign\F_GRCm38.81.P60Cerebellum_ALT.cluster.assign.RDS','\\iodine-cifs\broad_macosko\data\clusters\atlas_ica\F_GRCm38.81.P60Hippocampus\assign\F_GRCm38.81.P60Hippocampus.cluster.assign.RDS'};
%DropseqMeanAndVariancePath=fullfile(BeadseqCodePath,'DGEMeansAndVariances');

%NumPar=20; %number of threads
NumPar=1;

CropImage=0; %flag that asks you to crop output. THis is not used 
CropSuffix='_Cropped';

EnforceBaseBalance=1; 
BaseBalanceTolerance=0.05;

SaveData=1;
IlluminaReadThreshold=10;
MaximumBarcodesToAnalyze=160000;
NumLigations=20; %We assume that the missing ligation is Truseq-4 if you are doing 13.
NumBases=14; %This is the number of bases sequenced
BarcodeSequence=[1,2,3,4,0,5,0,6,0,7,8,9,10,11,0,12,0,13,0,14];
%BarcodeSequence=[1,2,3,4,0,5,0,6,7,8,9,10,0,11,0,12,0,13]; %this determines how the numerical barcode is built from the ligations. Basically, the constant bases should be 0, and all other bases should be consecutive
%tmpfolder='D:\pucktmp\';

ImageSize=6030; %The image registration will output images that are not all the exact same size, because of stitching. So find_roi_stack_fun crops the images a bit. We can choose this value to be something like 0.95 * the size of the images. So e.g. for 3x3 it should be 0.95*2048*(2*(2/3)+1) = 4500. For 7x7 it is 0.95*2048*(1+6*2/3)=9728. I used 10400 previously though.
XCorrBounds=[2800,3200,2800,3200]; %This is the ROI used in channel registration
RegisterColorChannels=1;
BeadZeroThreshold=1;
PixelCutoffRegistration=400;
PixelCutoffBasecalling=300;
DropBases=1;
BeadSizeCutoff=30;

%The illumina barcodes are assumed to be 13 bases long, consisting of the
%first 6 J bases and then the last 7 J bases. If the barcodes you are using
%are different, you have to change it in the bs2cs calls.

%And ligation sequnece:
%Tru_L1, Tru_L2, Tru-1_L1, Tru-1_L2, Tru-2_L1, Tru-2_L2, Tru-3_L1, Tru-3_L2, Tru-4_L1, Tru-4_L2
%The numbers in these variables are the ***second base being interrogated***
%Equivalently, it is the number associated with that ligation in figure 4a of the Fisseq nature protocols paper
PrimerNLigationSequence = [2, 7, 1, 6, 5, 4, 3]; %good for 14 ligations
%PrimerNLigationSequence = [2, 7, 1, 6, 5, 4]; %good for 13 ligations
%UP_L1, UP_L2, UP-1_L1, UP-1_L2, UP-2_L1, UP-2_L2, UP-3_L1, UP-3_L2, UP-4_L1, UP-4_L2
PrimerUPLigationSequence=[2, 7, 1, 6, 5, 4,3];

%Note that bs2cs will output the colors corresponding to each ligation in this ligation sequence, in the order specified here,
%So if you were only to do 6 ligations off primer N instead of 7, you would only need to
%remove the final element in the PrimerNLigationSequence

InverseLigationSequence=[3,1,7,6,5,4,2,10,8,14,13,12,11,9]; %Good for both 13 and 14 ligations.
%Before exporting, "N"s are added to the end of the ligation
%sequence, and must then be inserted into the correct location to stand in
%for the missing ligations. This code puts the N at position 7
%WhichLigationsAreMissing=[1,2,3,4,5,6,14,7,8,9,10,11,12,13];
WhichLigationsAreMissing=[1,2,3,4,5,6,7,8,9,10,11,12,13,14];

%Mask
Mask = logical([0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
%ogical mask for which ligation bases are bad and should not be used during
%mapping

%We need to determine if the device running the pipeline must be linked to
%the google doc for updating basecalling metrics into the puck log. Note
%that this must be set up for each unique device that runs the pipeline.
%The user will be prompted to manually enter a given security code into
%google.com/devices. After this is done once for a given device, subsequent
%pipeline runs can update automatically without any user action. Of note,
%privacy setting on puck log document must be set to "anyone with link can
%edit".

% reply = input('Is this the first time you will be running the Slide-seq pipeline on this device since 3/25/19 (reply: y/n)?','s');
% if strcmp(reply,'y')
%   %RunOnce(clientID,client_secret)
%   RunOnce('248080687373-vc0mevs6chk2rgpm44uva797tki6upf4.apps.googleusercontent.com', '6uYApyl6mlkwIk3cNH_qUVl')
% end


%% Load parameters from manifest file
display('Load parameters from manifest file');

% Input/select manifest file such as 'D:\Jilong\out\manifest.v2.txt'
% prompt = 'Please input the full path of your manifest file:\n';
% manifest = input(prompt);
%[file,path] = uigetfile('*.txt', 'Select manifest file');
%manifest = [path,file];

% Check if manifest file exists
if ~exist(manifest,'file')
    error('Manifest file not found');
end

% Load manifest content into a variable set
fid = fopen(manifest, 'rt');
fcon = textscan(fid,  '%s%s', 'Delimiter', '=');
fclose(fid);
indat = horzcat(fcon{:});
[m,n] = size(indat);
params = [];
for i=1:m
   params.(indat{i,1}) = indat{i,2}; 
end

% Retrieve variable values based on variable names
%FolderWithProcessedTiffs = params.('FolderWithProcessedTiffs');
%FolderWithProcessedTiffs = "/";
OutputFolderRoot = params.('OutputFolderRoot');
IlluminaFolderRoot = params.('IlluminaFolderRoot');
PucksToAnalyze = params.('PucksToAnalyze');
LibraryFolder = params.('LibraryFolder');
Reference = params.('Reference');
if isfield(params,'BeadType')
    BeadType = params.('BeadType');
end

PuckType="notSLAC";
if isfield(params,'PuckType')
    PuckType = params.('PuckType');
end

PuckName = params.('PuckName'); % Puck_190327
PuckNameStr = textscan(PuckName,'%s','Delimiter',',');
PuckNameStr = horzcat(PuckNameStr{:});
[m1,n1]=size(PuckNameStr);

PucksMatchedStr = params.('PucksMatched');
PucksMatchedStr = textscan(PucksMatchedStr,'%f','Delimiter',',');
PucksMatchedStr = horzcat(PucksMatchedStr{:});
[m,n]=size(PucksMatchedStr);
PucksMatched=string(m);
for i=1:m
	if m1==1
		PucksMatched(i)=[PuckName,'_',pad(num2str(PucksMatchedStr(i)),2,'left','0')];
	else
		PucksMatched(i)=[PuckNameStr{i},'_',pad(num2str(PucksMatchedStr(i)),2,'left','0')];
	end
end

% Convert string to cell array and vector
PucksToAnalyze = textscan(PucksToAnalyze,'%f','Delimiter',',');
PucksToAnalyze = horzcat(PucksToAnalyze{:});

% Retrieve Mask
if isfield(params,'Mask')
    Mask = textscan(params.('Mask'),'%f','Delimiter',',');
    Mask = transpose(horzcat(Mask{:}));
    Mask = logical(Mask);
end

if isfield(params,'IlluminaReadThreshold')
    IlluminaReadThreshold = textscan(params.('IlluminaReadThreshold'),'%f','Delimiter',',');
    IlluminaReadThreshold = transpose(horzcat(IlluminaReadThreshold{:}));
    IlluminaReadThreshold = double(IlluminaReadThreshold);
end

Monobase = 0;
if isfield(params,'Monobase')
    Monobase = textscan(params.('Monobase'),'%f','Delimiter',',');
    Monobase = transpose(horzcat(Monobase{:}));
    Monobase = logical(Monobase);
end

if Monobase == 1
	PrimerNLigationSequence = [1,2,3,4,5,6,7]; %good for 14 ligations
	%PrimerNLigationSequence = [2, 7, 1, 6, 5, 4]; %good for 13 ligations
	%UP_L1, UP_L2, UP-1_L1, UP-1_L2, UP-2_L1, UP-2_L2, UP-3_L1, UP-3_L2, UP-4_L1, UP-4_L2
	PrimerUPLigationSequence=[1,2,3,4,5,6,7];
	InverseLigationSequence=[1,2,3,4,5,6,7,8,9,10,11,12,13,14]; %Good for both 13 and 14 ligations.
end

% Create PuckNames from PuckName and PuckSToAnalyze
% Note that the order in PuckNames should match the order in the .nd2 file.
[m,n]=size(PucksToAnalyze);
PuckNames=string(m);
for i=1:m
	if m1==1
		PuckNames(i)=[PuckName,'_',pad(num2str(PucksToAnalyze(i)),2,'left','0')];
	else
		PuckNames(i)=[PuckNameStr{i},'_',pad(num2str(PucksToAnalyze(i)),2,'left','0')];
	end
end


%% Match Illumina barcodes - to rerun the bead mapping only, start here

display('Matching Illumina Barcodes')
%The illumina DGE should have the puck date and number, i.e. 180728_3,
%followed by an underscore, so 180728_3_
%The underscore is to distinguish 180728_2_ from 180728_23, for example

% Get folder info
OutputFoldersFinished='';
OutputFolders={};
for puck=1:length(PuckNames)
	if (PuckNames(puck)~=mypuck)
        continue
    end
    %ProcessedImageFolders{puck}=[FolderWithProcessedTiffs,PuckNames{puck},'/'];
    %OutputFolders{puck}=[OutputFolderRoot,PuckNames{puck},'/'];
    OutputFolders{puck}=[OutputFolderRoot,'/'];
    %mkdir([FolderWithProcessedTiffs,PuckNames{puck}]);
    %mkdir(OutputFolders{puck});
	OutputFoldersFinished=[OutputFolders{puck},'finished','/'];
end
if exist(OutputFoldersFinished,'dir')
	rmdir(OutputFoldersFinished);
end


MappingOutputFolders={};
for puck=1:length(PuckNames)
	if (PuckNames(puck)~=mypuck)
        continue
    end
	if Monobase == 1
        MappingOutputFolder=MapIlluminaToPuck_MaskBs(PuckNames{puck},PucksMatched(puck),BeadType,'PrimerNLigationSequence',PrimerNLigationSequence,'PrimerUPLigationSequence',PrimerUPLigationSequence,'InverseLigationSequence',InverseLigationSequence,'WhichLigationsAreMissing',WhichLigationsAreMissing,'IlluminaReadThreshold',IlluminaReadThreshold,'MaximumBarcodesToAnalyze',MaximumBarcodesToAnalyze,'OutputFolder',OutputFolders{puck},'IlluminaRootFolder',IlluminaFolderRoot,'NumBases',NumBases,'NumLigations',NumLigations,'BarcodeSequence',BarcodeSequence,'ImageSize',ImageSize,'NumPar',NumPar,"Mask",Mask,"LibraryFolder",LibraryFolder,"Reference",Reference,"PuckType",PuckType,"scriptfolder",scriptfolder);
    else
        MappingOutputFolder=MapIlluminaToPuck_Mask(PuckNames{puck},PucksMatched(puck),BeadType,'PrimerNLigationSequence',PrimerNLigationSequence,'PrimerUPLigationSequence',PrimerUPLigationSequence,'InverseLigationSequence',InverseLigationSequence,'WhichLigationsAreMissing',WhichLigationsAreMissing,'IlluminaReadThreshold',IlluminaReadThreshold,'MaximumBarcodesToAnalyze',MaximumBarcodesToAnalyze,'OutputFolder',OutputFolders{puck},'IlluminaRootFolder',IlluminaFolderRoot,'NumBases',NumBases,'NumLigations',NumLigations,'BarcodeSequence',BarcodeSequence,'ImageSize',ImageSize,'NumPar',NumPar,"Mask",Mask,"LibraryFolder",LibraryFolder,"Reference",Reference,"PuckType",PuckType,"scriptfolder",scriptfolder);    
    end
	if (MappingOutputFolder=="NA")
        continue
    end
    
    % Copy manifest file to mapping output folder
    copyfile(manifest,[MappingOutputFolder,'/']);
    
    MappingOutputFolders{puck}=MappingOutputFolder;
    %Record parameters into a file, "IlluminaParams.txt"
    fileid=fopen([MappingOutputFolder,'IlluminaParams.txt'],'a');
    fprintf(fileid,['\n\nThe parameters for this Bead-To-Illumina matching run were:\n',...
    '\nBead Type: ',char(BeadType),...
    '\nPrimerNLigationSequence: ',num2str(PrimerNLigationSequence),...
    '\nPrimerUPLigationSequence: ',num2str(PrimerUPLigationSequence),...
    '\nInverseLigationSequence: ',num2str(InverseLigationSequence),...
    '\nWhichLigationsAreMissing: ',num2str(WhichLigationsAreMissing),...
    '\nIlluminaReadThreshold: ',num2str(IlluminaReadThreshold),...
    '\nMaximumBarcodesToAnalyze: ',num2str(MaximumBarcodesToAnalyze),...
    '\nNumBases: ',num2str(NumBases),...
    '\nNumLigations: ',num2str(NumLigations),...
    '\nBarcodeSequence: ',num2str(BarcodeSequence)]);
    fclose(fileid);
%    ConvertMatToR(OutputFolders{puck});
%    tmp=strsplit(MappingOutputFolders{puck},'\');
%    mappingstarttimereadable=tmp(end-1);
%    mappingstarttimereadable=mappingstarttimereadable{1};
%    mkdir(['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',mappingstarttimereadable]);
%    copyfile(MappingOutputFolders{puck},['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',mappingstarttimereadable]);
end

close all

if CropImage
    disp('Cropping images')
    for puck=1:length(PuckNames)
		if (PuckNames(puck)~=mypuck)
			continue
		end
        disp(['Loading files to crop puck',num2str(puck)])
        try
            BeadMappingFile=FindMostRecentMapping(OutputFolders{puck});
        catch
            continue
        end
        MappingOutputFolder=fullfile(OutputFolders{puck},BeadMappingFile);
        if ~exist(fullfile(MappingOutputFolder,'IlluminaParams.txt')) %Puck hasn't been mapped
            continue
        end
        load(fullfile(MappingOutputFolder,'BijectiveMapping.mat'))
        image=PlotGeneFromName('Malat1',GeneNames,UniqueMappedDGE,UniqueMappedBeads,'Overlay',1,'PlotStyle','Default');
        disp(['Ready to crop puck',num2str(puck)])
        h = imfreehand; %draw something 
        M = h.createMask();
        
        disp('Cropping complete')
        Locs=[UniqueMappedBeads.Locations];
%        GoodBeads=Locs(1,:) > croprectangle(1) & Locs(1,:) < (croprectangle(1)+croprectangle(3)) & Locs(2,:) > croprectangle(2) & Locs(2,:)<croprectangle(2)+croprectangle(4);
        GoodBeads=M(sub2ind(size(image),round(Locs(2,:)),round(Locs(1,:))))==1;
        UniqueMappedBeads=UniqueMappedBeads(GoodBeads);
        UniqueMappedDGE=UniqueMappedDGE(:,GoodBeads);
        UniqueMappedIlluminaBarcodes=UniqueMappedIlluminaBarcodes(GoodBeads);
%        movefile(fullfile(MappingOutputFolder,'BeadLocationsForR.csv'),fullfile(MappingOutputFolder,'BeadLocationsForRUncropped.csv'));
%        movefile(fullfile(MappingOutputFolder,'MappedDGEForR.csv'),fullfile(MappingOutputFolder,'MappedDGEForRUncropped.csv'));    
%        movefile(fullfile(MappingOutputFolder,'BijectiveMapping.mat'),fullfile(MappingOutputFolder,'BijectiveMappingUncropped.mat'));    
        save(fullfile(MappingOutputFolder,['BijectiveMapping',CropSuffix,'.mat']),'UniqueMappedBeads','UniqueMappedDGE','UniqueMappedIlluminaBarcodes','GeneNames','-v7.3');
        ConvertMatToR(OutputFolders{puck},'CropSuffix',CropSuffix);    
    end
end


mkdir(OutputFoldersFinished);


% for puck=1:length(PuckNames)
    % if ~RunAnalogizer(puck)
        % continue
    % end
    % AtlasReference=RunAnalogizer(puck);
    % if AnalogizerType=="Analogizer"
        % DropseqDGEPath=DropseqDGEPaths{AtlasReference};
        % DropseqClusterPath=DropseqClusterPaths{AtlasReference};
    % else
        % DropseqDGEPath='Null';
        % DropseqClusterPath='Null';
    % end
    % display('Running Analogizer')
    % %This will run the analogizer on the most recent Illumina Barcoding run.
    % try
        % BeadMappingFile=FindMostRecentMapping(OutputFolders{puck});
    % catch
        % continue
    % end

    % MappingOutputFolder=fullfile(OutputFolders{puck},BeadMappingFile);
    % if ~exist(fullfile(MappingOutputFolder,'IlluminaParams.txt')) %Puck hasn't been mapped
        % continue
    % end

    % if exist(fullfile(MappingOutputFolder,'AnalogizerParams.txt'))
        % c=clock;
        % analogizerstarttimereadable=[num2str(c(2)),'-',num2str(c(3)),'_',pad(num2str(c(4)),2,'left','0'),pad(num2str(c(5)),2,'left','0')];

        % %If prior analogizer data exists, we save it and its parameter file
        % %into a new labeled folder.
        % mkdir(fullfile(MappingOutputFolder,['DataFromAnalogizer_MovedOn_',analogizerstarttimereadable]))
        % try
            % try
                % movefile(fullfile(MappingOutputFolder,'Cluster*'),fullfile(MappingOutputFolder,['DataFromAnalogizer_MovedOn_',analogizerstarttimereadable]))            
            % catch
            % end
            % try
                % movefile(fullfile(MappingOutputFolder,'Analogizer*'),fullfile(MappingOutputFolder,['DataFromAnalogizer_MovedOn_',analogizerstarttimereadable]))            
            % catch
            % end
            % %This is a failsafe
            % try
                % copyfile(fullfile(MappingOutputFolder,'Analogizer*'),fullfile(MappingOutputFolder,['DataFromAnalogizer_MovedOn_',analogizerstarttimereadable]))
            % catch
            % end
            % try
                % copyfile(fullfile(MappingOutputFolder,'Cluster*'),fullfile(MappingOutputFolder,['DataFromAnalogizer_MovedOn_',analogizerstarttimereadable]))
            % catch
            % end

            % delete(fullfile(MappingOutputFolder,'Analogizer*'))
            % delete(fullfile(MappingOutputFolder,'AnalogizerParams.txt'))
            % delete(fullfile(MappingOutputFolder,'Cluster*'))
        % catch

        % end
    % end
    % thisbeadcutoff=AnalogizerBeadCutoff;
    
    % while ~exist(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'))
        % commandfile=fopen('C:\Analogizer.cmd','w');
        % %fwrite(commandfile,['cd "C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code" & dir']);
        % if AnalogizerType=="Analogizer"
            % fwrite(commandfile,['cd "C:\Users\sgr\Dropbox (MIT)\Project - SlideSeq\BeadSeq Code" & "C:\Program Files\R\R-3.4.3\bin\Rscript.exe" AnalogizerScript.R ',PuckNames{puck},'\',BeadMappingFile,' ',DropseqDGEPath,' ',DropseqClusterPath,' ',num2str(thisbeadcutoff)]);
        % elseif AnalogizerType=="NMFReg"
            % %Bob says:
% %           -da = data path atlas
% %           -dp = data path puck
% %           -t = tissue type (only cerebellum and hippo now
% %           - c = UMI cutoff
% %           - dge = DGE name minus file extension
% %           - bl = Bead locations minus file extension
            % switch RunAnalogizer(puck)
                % case 1
                    % tissuetype='cerebellum';
                % case 2
                    % tissuetype='hippocampus';
                % case 3
                    % tissuetype='frontalcortex';
                % case 4
                    % tissuetype='posteriorcortex';
                % otherwise
                    % assert(1==0)
            % end
            % fwrite(commandfile,['C: & cd "',BeadseqCodePath,'" & "',PythonPath,'" autoNMFreg_windows.py -da \\iodine-cifs\broad_macosko\data\NMFreg\data -dp "',fullfile(OutputFolders{puck},BeadMappingFile),'" -t ',tissuetype,' -c ',num2str(thisbeadcutoff),' -dge MappedDGEForR -bl BeadLocationsForR']);            
        % end
        % fclose(commandfile);
        % !C:/Analogizer
        % if ~exist(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'))
            % disp('Analogizer Failed. Increasing bead cutoff and continuing.')
            % thisbeadcutoff=thisbeadcutoff+5;
        % end
        % if AnalogizerType=="NMFReg" && exist(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv')) && ~exist(fullfile(MappingOutputFolder,'AnalogizerClusterAssignmentsOriginal.csv'))
            % %NMFReg outputs the cluster assignments in the wrong format
            % opts=detectImportOptions(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'));
            % ClusterAssignments=readtable(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'),opts,'ReadVariableNames',true);
            % newtable=table(ClusterAssignments.barcode,ClusterAssignments.atlas_cluster,'VariableNames',{'Var1','x'});
            % movefile(fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'),fullfile(MappingOutputFolder,'AnalogizerClusterAssignmentsOriginal.csv'))
            % writetable(newtable,fullfile(MappingOutputFolder,'AnalogizerClusterAssignments.csv'));
        % end
    % end
    % switch RunAnalogizer(puck)
        % case 1
            % copyfile(fullfile(BeadseqCodePath,'\DGEMeansAndVariances\CerebellumVarianceByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerVarianceByDropseqCluster.csv'));
            % copyfile(fullfile(BeadseqCodePath,'\DGEMeansAndVariances\CerebellumExpressionByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerExpressionByDropseqCluster.csv'));
        % case 2
            % copyfile(fullfile(BeadseqCodePath,'\DGEMeansAndVariances\HippocampusVarianceByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerVarianceByDropseqCluster.csv'));
            % copyfile(fullfile(BeadseqCodePath,'\DGEMeansAndVariances\HippocampusExpressionByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerExpressionByDropseqCluster.csv'));
        % case 3
            % copyfile(fullfile(BeadseqCodePath,'\DGEMeansAndVariances\FrontalCortexVarianceByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerVarianceByDropseqCluster.csv'));
            % copyfile(fullfile(BeadseqCodePath,'\DGEMeansAndVariances\FrontalCortexExpressionByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerExpressionByDropseqCluster.csv'));
        % case 4
            % copyfile(fullfile(BeadseqCodePath,'\DGEMeansAndVariances\PosteriorCortexVarianceByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerVarianceByDropseqCluster.csv'));
            % copyfile(fullfile(BeadseqCodePath,'\DGEMeansAndVariances\PosteriorCortexExpressionByDropseqCluster.csv'),fullfile(MappingOutputFolder,'AnalogizerExpressionByDropseqCluster.csv'));    
    % end
    % %Record parameters into a file, "AnalogizerParams.txt"
    % fileid=fopen([MappingOutputFolder,'\AnalogizerParams.txt'],'a');
    % fprintf(fileid,['\n\nThe parameters for this Bead-To-Illumina matching run were:\n',...
    % '\nAnalogizer Type: ',char(AnalogizerType),...
    % '\nAnalogizer Bead Cutoff: ',num2str(thisbeadcutoff),...
    % '\nDropseq DGE Path: ',DropseqDGEPath,...
    % '\nDropseq Cluster Path: ',DropseqClusterPath,...
    % ]);
    % fclose(fileid);
    % if RunSignificanceAnalysis
        % for k=1:NumClusters(AtlasReference)
            % PermutationTestByCluster(OutputFolders{puck},k)
        % end
    % end
% end

%for puck=1:length(PuckNames)  
%    tmp=strsplit(MappingOutputFolders{puck},'\');
%    mappingstarttimereadable=tmp(end-1);
%    mappingstarttimereadable=mappingstarttimereadable{1};
%    mkdir(['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',mappingstarttimereadable]);
%    copyfile(MappingOutputFolders{puck},['\\iodine-cifs\broad_macosko\data\pucks\',PuckNames{puck},'_',mappingstarttimereadable]);
%end

%% Evaluate

%We can now use the Hough transform, and ask what percentage of Hough beads
%also have barcodes, as a way of analyzing what percentage of the beads
%were called. We call the Hough transform and find the centroids, and then
%ask about the fraction of Hough beads with centroids within a radius of a sequenced
%barcode ;

%We want to produce a plot with the fraction of illumina barcodes that are diporect
%matches with a surface bead; one base away from a surface bead; two bases,
%etc., and also for the negative control
