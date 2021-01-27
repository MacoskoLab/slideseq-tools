function WriteBijectiveMapping(BeadFile, DgeFile, IlluminaBarcodeFile, GeneNameFile, LocationFile, OutputFolder)

c=clock;
mappingstarttimereadable=[num2str(c(1)),'_',num2str(c(2)),'-',num2str(c(3)),'_',pad(num2str(c(4)),2,'left','0'),pad(num2str(c(5)),2,'left','0')];
mkdir([OutputFolder,'/','BeadMapping_',mappingstarttimereadable]);
MappingOutputFolder=[OutputFolder,'/','BeadMapping_',mappingstarttimereadable,'/'];

BeadBarcodes = fileread(BeadFile);
BeadBarcodes = strsplit(BeadBarcodes);
BeadBarcodes(end)=[];

UniqueMappedIlluminaBarcodes = fileread(IlluminaBarcodeFile);
UniqueMappedIlluminaBarcodes = strsplit(UniqueMappedIlluminaBarcodes);
UniqueMappedIlluminaBarcodes(end)=[];

GeneNames = fileread(GeneNameFile);
GeneNames = strsplit(GeneNames);
GeneNames(end)=[];

BeadLocations = dlmread(LocationFile);
BeadLocations(:,1) = [];
BeadLocations = transpose(BeadLocations);
UniqueMappedBeads=struct('Barcodes',num2cell(BeadBarcodes),'Locations',num2cell(BeadLocations,1));

UniqueMappedDGE = dlmread(DgeFile);

save([MappingOutputFolder,'BijectiveMapping.mat'],'UniqueMappedBeads','UniqueMappedDGE','UniqueMappedIlluminaBarcodes','GeneNames','-v7.3');
