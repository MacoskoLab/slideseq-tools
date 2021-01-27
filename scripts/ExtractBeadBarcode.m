function ExtractBeadBarcode(PuckCallerFolder)

NumBases=14; %This is the number of bases sequenced
WhichLigationsAreMissing=[1,2,3,4,5,6,7,8,9,10,11,12,13,14];

load([PuckCallerFolder,'AnalysisOutputs-selected'],'Bead');
    
BeadBarcodes = [Bead.Barcodes];
BeadLocations = {Bead.Locations};
[UniqueBeadBarcodes,BBFirstRef,BBOccCounts]=unique(BeadBarcodes, 'stable');
[~, b] = ismember(BeadBarcodes, BeadBarcodes(BBFirstRef(accumarray(BBOccCounts,1)>1)));
UniqueBeadBarcodes2 = BeadBarcodes(b<1);
UniqueBeadLocations = BeadLocations(b<1);
BaseBalanceBarcodes=[UniqueBeadBarcodes2];
BaseBalanceBase5Barcodes=cellfun(@(x) reverse(string(x)),{dec2base(BaseBalanceBarcodes,5,NumBases)},'UniformOutput',false);
BaseBalanceBase5Barcodes=BaseBalanceBase5Barcodes{1};
UniqueBeadBarcodesForExport=char(replace(BaseBalanceBase5Barcodes,{'0','1','2','3','4'},{'N','T','G','C','A'}));
if NumBases<14 %This is to deal with the InverseLigationSequence -- the export barcodes have to be 14 bases long
	UniqueBeadBarcodesForExport(:,NumBases+1:14)='N';
	UniqueBeadBarcodesForExport=UniqueBeadBarcodesForExport(:,WhichLigationsAreMissing);
end
file=fullfile(PuckCallerFolder,'BeadBarcodes.txt');
dlmwrite(file,[UniqueBeadBarcodesForExport]);
file=fullfile(PuckCallerFolder,'BeadLocations.txt');
dlmwrite(file,[UniqueBeadLocations]);
