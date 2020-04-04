disp(['Running LungCancerSaveData.m']) ;
%  Does data organization and saving of 
%  Early TCGA Data with RNAseq read depth curves for 4 genes. 
%  Data from Matt Wilkerson, Oct. 2011. 
%
%  Inputs in directory OriginalData:
%  counts.csv:  Main file of read counts. Top row is case ID, 1st column is chromosome location, remaining columns are read counts
%  exonsMarron.csv:  Edited version, to fix errors in original file exons.csv  (by deleting the lines with the 1st column = 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 54.  Look at P10, Case Number 38)
%
%  This script is a modification of Lung Cancer 2011 
%          (posted on STOR 881 webpage 2017 & 2019),
%      which mostly improves it by considering all 4 genes
%


datfilename = '../DerivedData/LungCancer2011' ;


%  First Read in Exon Info
%
filestr = '../OriginalData/exonsMarron.csv' ;
[numeric,txt]=xlsread(filestr) ;
vchromall = numeric(:,2) ;
vexleft = numeric(:,3) ;
vexright = numeric(:,4) ;

ng = 4 ;    %  Number of genes in input spreadsheet

vchrom = unique(vchromall) ;

disp('  Check chromosome numbers:    10 19 3 9') ;
vchrom'

GeneNamesS = {'PIK3CA'; 'CDKN2A'; 'P10'; 'STK11'; } 
    %  Structure with Gene Names
    %  Manually copied from spreadsheet
    %  Put in same order as in vchrom
    %      Note CDK2A was used in the file, but should really be CDKN2A

disp('  Check gene names against entries in file:') ;
for ig = 1:ng ;
  disp(['    For Chromosome ' num2str(vchrom(ig)) ...
                    '    Gene Name is ' GeneNamesS{ig}]) ;
end ;



disp(' ') ;

%  Now create structure of exon lefts and rights
%
vexleftS = {} ;
vexrightS = {} ;
for ig = 1:ng ;
  flag = (vchromall == vchrom(ig)) ;
  vexleftS{ig} = unique(vexleft(flag)) ;
  vexrightS{ig} = unique(vexright(flag)) ;
end ;

disp('  Check 1st exon left is 178866311 by showing difference is 0') ;
abs(vexleftS{1}(1) - 178866311)

disp('  Check last exon left is 1227592 by showing difference is 0') ;
abs(vexleftS{end}(end) - 1227592)

disp('  Check 1st exon right is 178866391 by showing difference is 0') ;
abs(vexrightS{1}(1) - 178866391)

disp('  Check last exon right is 1228434 by showing difference is 0') ;
abs(vexrightS{end}(end) - 1228434)



%  Now Read in Main Data
%
filestr = '../OriginalData/counts.csv' ;
[numeric,txt]=xlsread(filestr) ;
vbpnall = numeric(:,1) ;
    %  all base pair numbers in input file
mctsall = numeric(:,2:end) ;
    %  matrix of all counts
CaseNamesS = txt ;
    %  Structure of Case Names

numeric = [] ;
txt = [] ;
    %  to save space


%  Check Inputs
%
disp('  ') ;

disp('Check 1st Base Pair is 1205795 by showing difference is 0') ;
abs(vbpnall(1) - 1205795)

disp('Check last Base Pair is 89728533 by showing difference is 0') ;
abs(vbpnall(end) - 89728533)

disp('Check Upper Left Count is 0') ;
mctsall(1,1)

disp('Check Upper Right Count is 0') ;
mctsall(1,end)

disp('Check Lower Left Count is 29') ;
mctsall(end,1)

disp('Check Lower Right Count is 7') ;
mctsall(end,end)

disp('Check 1st Case Name is preview.20110919/110107_UNC9-SN296_0118_B815C0ABXX7TCGA-66-2756-11A-01R-A') ;
CaseNamesS{1}

disp('Check last Case Name is preview.20110919/110630_UNC11-SN627_0112_AD0CVJABXX8TCGA-22-5482-01A-01R-1635-07') ;
CaseNamesS{end}


%  Put Data into Structure
%
DataS = {} ;
for ig = 1:ng ;

  vexleft = vexleftS{ig} ;
  vexright = vexrightS{ig} ;

  %  Check vectors have the same length
  %
  if length(vexleft) ~= length(vexright) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Error from VisualizeNextGen2011:   !!!') ;
    disp('!!!   Uneven length of Exon              !!!') ;
    disp('!!!   Boundary Vectors                   !!!') ;
    disp('!!!   Terminating Execution              !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    return ;
  end ;

  %  Check each left exon end is larger than each right exon end
  %
  if sum(vexright <= vexleft) > 0 ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Error from VisualizeNextGen2011:   !!!') ;
    disp('!!!   A left exon boundary is smaller    !!!') ;
    disp('!!!   than the right boundar             !!!') ;
    disp('!!!   Terminating Execution              !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    return ;
  end ;


  vbpn = [] ;
  mcts = [] ;
  for il = 1:length(vexleft) ;
    flag = (vexleft(il) <= vbpnall) & (vbpnall <= vexright(il)) ;
    vbpn = [vbpn; vbpnall(flag)] ;
    mcts = [mcts; mctsall(flag,:)] ;
  end ;


  DataS{ig,1} = vbpn ;
  DataS{ig,2} = mcts ;

end ;

vbpnall = [] ;
mctsall = [] ;
    %  to save space


%  Recreate  first page of Matt's file: stk11.summary.pdf
%
figure(1) ;
clf ;
ig = 3 ;    %  for Chromosome 10 = P10
plot(DataS{ig,2}(:,1),'k-') ;
title(CaseNamesS{1}) ;
xlabel('exonic nt number, not genomic position') ;
ylabel('RNA read depth') ;
vaxh = axisSM([1:size(DataS{ig,2},1)]) ;
vaxv = [0; 610] ;
text(vaxh(1) + 0.1 * (vaxh(2) - vaxh(1)), ...
     vaxv(1) + 0.9 * (vaxv(2) - vaxv(1)), ...
     ['Chromosome ' num2str(vchrom(ig)) '    Gene = ' GeneNamesS{ig}]) ;
axis([vaxh(1) vaxh(2) vaxv(1) vaxv(2)]) ;


%  Recreate  last page of Matt's file: stk11.summary.pdf
%
figure(2) ;
clf ;
ig = 2 ;    %  for Chromosome 9 = CDKN2A
plot(DataS{ig,2}(:,end),'k-') ;
title(CaseNamesS{end}) ;
xlabel('exonic nt number, not genomic position') ;
ylabel('RNA read depth') ;
vaxh = axisSM([1:size(DataS{ig,2},1)]) ;
vaxv = [0; 2010] ;
text(vaxh(1) + 0.1 * (vaxh(2) - vaxh(1)), ...
     vaxv(1) + 0.9 * (vaxv(2) - vaxv(1)), ...
     ['Chromosome ' num2str(vchrom(ig)) '    Gene = ' GeneNamesS{ig}]) ;
axis([vaxh(1) vaxh(2) vaxv(1) vaxv(2)]) ;



%  Save Data as .mat File
%
save(datfilename,'DataS','CaseNamesS','vchrom','GeneNamesS') ;
    %  Save variables:
    %      DataS
    %      CaseNamesS
    %      vchrom
    %      GeneNamesS
disp(['  File ' datfilename '.mat finished writing']) ; 
disp(' ') ;


%  Save Data as .xls Files, one for each gene,
%  to be turned into .csv files in Excel.
%
for ig = 1:4 ;

  disp(' ') ;
  dfname = [datfilename GeneNamesS{ig}] ;

  status0 = xlswrite(dfname,{'Genomic Coords'},1,'A1') ;
  if status0 ;
    disp(['  Wrote "Genomic Coordinates" to ' dfname '.xls']) ;
  else ;
    disp(['  Column headers write to ' dfname '.xls failed']) ;
  end ;

  status1 = xlswrite(dfname,CaseNamesS,1,'B1') ;
  if status1 ;
    disp(['  Wrote Column headers to ' dfname '.xls']) ;
  else ;
    disp(['  Column headers write to ' dfname '.xls failed']) ;
  end ;

  status2 = xlswrite(dfname,DataS{ig,1},1,'A2') ;
  if status2 ;
    disp(['  Wrote Row headers to ' dfname '.xls']) ;
  else ;
    disp(['  Row headers write to ' dfname '.xls failed']) ;
  end ;

  status3 = xlswrite(dfname,DataS{ig,2},1,'B2') ;
  if status3 ;
    disp(['  Wrote Main Data to ' dfname '.xls']) ;
  else ;
    disp(['  Main data write to ' dfname '.xls failed']) ;
  end ;

  if (status1 & status2 & status3) ;
    disp(['Should now save ' dfname '.xls to a .csv file in Excel']) ;
  end ;

end ;




