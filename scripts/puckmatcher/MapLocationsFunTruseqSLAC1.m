function Illumina=MapLocationsFunTruseqSLAC1(OutputDirectory,Bead,IlluminaBarcodes,PrimerNLigationSequence,NumBases,SaveData,PuckName)
%The goal of this function is:
%1) identify bead barcodes with illumina barcodes
%2) QC that identification: make the histogram of hamming distances for the
%test case and the negative control.
%3) Output a reordered form of BeadLocations, where BeadLocation(k) is the
%location for the kth bead in the DGE.

%This version of MapLocationsFun assumes the following bead structure:
%T JJJJJJJJJJJJ
%Note that we only have one primer here.

%NOTE: We convert the Illumina barcodes into color space, because of the
%strange error properties of the dibase encoding in SOLiD. Thus, we do not
%have to reconstruct the base-space sequence associated with the SOLiD
%sequencing explicitly. So, it doesn't matter if we don't have enough
%ligations to reconstruct the base-space sequence -- the comparison here
%should still work.


%% Base space to color space conversion
%This is a super shitty thing we have to do because of the way we're doing
%the ligations

%for UP, we do primer UP ligation 1, then primer UP-1 ligation 1,
% Primer UP-1 ligation 2, primer UP-2 ligation 2, primer UP-3 ligation 2, and
% primer UP-4 ligation 2

Base5Barcodes=zeros(1,length(IlluminaBarcodes));
NegBarcodes=zeros(1,length(IlluminaBarcodes));
badflags=0; 
for seqnum=1:length(IlluminaBarcodes)
    seq=char(IlluminaBarcodes(seqnum));
    if length(seq)<12 || string(seq)=="NNNNNNNNNNNN" %this in particular messes us up, because then base5barcodes=0.
        Base5Barcodes(seqnum)=1;
        NegBarcodes(seqnum)=1;
        continue
    end
    [PrimerNcolorspace,badflagN]=bs2cs(seq(1:12),PrimerNLigationSequence,'T','');
    colorspacesequence=PrimerNcolorspace;
    if badflagN
        badflags=badflags+1;
    end
    randompermutation=randperm(NumBases);
    for m=1:length(colorspacesequence)% We just do 1:11 right now because puck 3 is missing the last ligation
        Base5Barcodes(seqnum)=Base5Barcodes(seqnum)+uint64(colorspacesequence(m))*5^(m-1);
        NegBarcodes(seqnum)=NegBarcodes(seqnum)+uint64(colorspacesequence(m))*5^(randompermutation(m)-1);
    end
end
%For Primer N, we do Primer N ligation 1, Primer N ligation 2, Primer N-1 ligation 1, Primer N-1 ligation 2, Primer N-2 ligation 2, Primer N-3 ligation 2 
save([OutputDirectory,'Base5Barcodes'],'Base5Barcodes');

%% Match Barcodes
%We now match the SOLiD barcodes with the IlluminaBarcodes. Note that this
%assumes that we are doing parallelization over 20 nodes.
    [UniqueBeadBarcodes,BBFirstRef,BBOccCounts]=unique([Bead.Barcodes]);
    UniqueBeadLocations={Bead(BBFirstRef).Locations};
    UniqueBeadPix={Bead(BBFirstRef).Pixels};
    %    [IdentifiedBarcodes,IA,IB]=intersect(BeadBarcodes, Base5Barcodes);
%    [negIdentifiedBarcodes,negIA,negIB]=intersect(BeadBarcodes, NegBarcodes);
    
    paddedBase5Barcodes=zeros(1,ceil(length(Base5Barcodes)/20)*20);
    paddedBase5Barcodes(1:length(Base5Barcodes))=Base5Barcodes;
    paddedNegBarcodes=zeros(1,ceil(length(NegBarcodes)/20)*20);
    paddedNegBarcodes(1:length(NegBarcodes))=NegBarcodes;
    reshapedBase5Barcodes=reshape(paddedBase5Barcodes,ceil(length(Base5Barcodes)/20),20);
    reshapedNegBarcodes=reshape(paddedNegBarcodes,ceil(length(NegBarcodes)/20),20);
    hammingdistancecell={};
    neghammingdistancecell={};
    IndexofBeadBarcodecell={};
    NumNearestBarcodescell={};
    delete(gcp('nocreate'))
    pool=parpool(20);
    
    parfor parnum=1:20 %this could maybe be accelerated by sorting the barcodes?
        localNegBarcodes=reshapedNegBarcodes(:,parnum);
        localBase5Barcodes=reshapedBase5Barcodes(:,parnum);
        localIndexofBeadBarcode=zeros(1,nnz(localBase5Barcodes));
        localNumNearestBarcodes=zeros(1,nnz(localBase5Barcodes)); %when we find the minimum hamming distance, this is the 
        localneghammingdistances=zeros(1,nnz(localNegBarcodes));
        localhammingdistances=zeros(1,nnz(localBase5Barcodes));
    for seq=1:nnz(localBase5Barcodes)
        if seq/1000==ceil(seq/1000)
            disp(['Worker ',num2str(parnum),' is on barcode ',num2str(seq)])
        end
        hammingdistancetmp=cellfun(@(x) dec2base(x,5,NumBases)==dec2base(localBase5Barcodes(seq),5,NumBases), {UniqueBeadBarcodes},'UniformOutput',false);
        tmpsum=sum(cell2mat(hammingdistancetmp),2);        
        [mval,ival]=max(tmpsum);
        localhammingdistances(seq)=NumBases-mval;
        localIndexofBeadBarcode(seq)=ival;
        localNumNearestBarcodes(seq)=sum(tmpsum==mval);
%        if mval==14 && localNumNearestBarcodes(seq)>1
%            display(['problem with barcode ',dec2base(localBase5Barcodes(seq),5,14)]);
%        end
        neghammingdistancetmp=cellfun(@(x) dec2base(x,5,NumBases)==dec2base(localNegBarcodes(seq),5,NumBases), {UniqueBeadBarcodes},'UniformOutput',false);
        localneghammingdistances(seq)=NumBases-max(sum(cell2mat(neghammingdistancetmp),2));
%        hammingdistances(seq)=min(cellfun(@(x) sum(dec2base(abs(x-Base5Barcodes(seq)),5,14)=='1'),{BeadBarcodes},'UniformOutput',false));
%        min=15;
%        for seq2=1:length(BeadBarcodes)
%            num=sum(dec2base(abs(BeadBarcodes(seq2)-Base5Barcodes(seq)),5,14)=='1');
%            if num<min
%                min=num;
%            end
%        end
    end
    hammingdistancecell{parnum}=localhammingdistances;
    neghammingdistancecell{parnum}=localneghammingdistances;
    IndexofBeadBarcodecell{parnum}=localIndexofBeadBarcode;
    NumNearestBarcodescell{parnum}=localNumNearestBarcodes;
    end
    delete(pool);
    numseqs=0;
    for ppp=1:20
        numseqs=numseqs+length(hammingdistancecell{ppp});
    end
    hammingdistances=zeros(1,numseqs);
    neghammingdistances=zeros(1,numseqs);
    IndexofBeadBarcode=zeros(1,numseqs);
    NumNearestBarcodes=zeros(1,numseqs);
    HammingDistanceIndex=1;
    for k=1:20
        hammingdistances(HammingDistanceIndex:(HammingDistanceIndex+length(hammingdistancecell{k})-1))=hammingdistancecell{k};
        neghammingdistances(HammingDistanceIndex:(HammingDistanceIndex+length(neghammingdistancecell{k})-1))=neghammingdistancecell{k};
        IndexofBeadBarcode(HammingDistanceIndex:(HammingDistanceIndex+length(hammingdistancecell{k})-1))=IndexofBeadBarcodecell{k};
        NumNearestBarcodes(HammingDistanceIndex:(HammingDistanceIndex+length(hammingdistancecell{k})-1))=NumNearestBarcodescell{k};
        HammingDistanceIndex=HammingDistanceIndex+length(hammingdistancecell{k});
    end
    
    Illumina=struct('Barcode',IlluminaBarcodes,'HammingDistance',num2cell(hammingdistances),'NumNearestBarcodes',num2cell(NumNearestBarcodes),'Base5Barcode',num2cell(Base5Barcodes),'IndexofBeadBarcode',num2cell(IndexofBeadBarcode),'SOLiDBarcode',num2cell(UniqueBeadBarcodes(IndexofBeadBarcode)),'MappedLocation',UniqueBeadLocations(IndexofBeadBarcode),'Pixels',UniqueBeadPix(IndexofBeadBarcode));
    
%    MappedLocations=UniqueBeadLocations(:,IndexofBeadBarcode);
    if SaveData
        save([OutputDirectory,'NegBarcodes.mat'],'NegBarcodes','neghammingdistances');
        save([OutputDirectory,'Illumina.mat'],'Illumina','-v7.3');
    end
%    save([OutputDirectory,'MappedLocations.mat'],'MappedLocations');
%    save([OutputDirectory,'hammingdistances.mat'],'hammingdistances');
    
    figure(5)
    clf
    histogram(hammingdistances,0:1:NumBases)
    hold on
    histogram(neghammingdistances,0:1:NumBases)
    title('Minimum Hamming Distance Between Beads and Illumina Barcodes or Negative Control')
    export_fig([OutputDirectory,'Report_',PuckName,'.pdf'],'-append');
