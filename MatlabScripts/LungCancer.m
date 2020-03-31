disp(['Running LungCancer.m']) ;
%  Does data organization and simple visualization of 
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

ipart = 1 ;
              %    0 - Read in Raw Data, and Save as .mat file
              %    1 - Raw Data Curves, original scale
              %    2 - Raw Data Curves, log10 scale
              %    3 - PCA scatterplot
              %    4 - PCA scatterplot, brushed
              %    5 - Raw Data Curves, log10 scale, brushed


datfilename = '../DerivedData/LungCancer2011' ;


if ipart == 0 ;    %  Read in Raw Data, and Save as .mat file

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

  GeneNamesS = {'PIK3CA'; 'CDK2A'; 'P10'; 'STK11'; } 
      %  Structure with Gene Names
      %  Manually copied from spreadsheet
      %  Put in same order as in vchrom
      %      Note CDK2A is used in the file, but should really be CDKN2A

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
  if ipart == 0 ;
    filestr = '../OriginalData/counts.csv' ;
  end ;
  [numeric,txt]=xlsread(filestr) ;
  vbpnall = numeric(:,1) ;
      %  all base pair numbers in input file
  mctsall = numeric(:,2:end) ;
      %  matrix of all counts
  CaseNamesS = txt(1,2:end) ;
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
  ig = 2 ;    %  for Chromosome 9 = CDK2A
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



else ;    %  Load data from previously saved .mat file

  %  Load data from .mat file
  %
  load(datfilename) ;
      %  Loads variables:
      %      DataS
      %      CaseNamesS
      %      vchrom
      %      GeneNamesS


  %  Set basics
  %
  ig = 2 ;    %  for Chromosome 9 = CDK2NA
  n = length(CaseNamesS) ;


  mdat = DataS{ig,2} ;
      %  raw data, for gene CDKN2A

  mdatl = log10(mdat + 1) ;


  %  Set up brushing colors
  %
  %  Based on initial PCA
  %
  paramstruct = struct('npc',4,...
                       'iscreenwrite',1,...
                       'viout',[0 0 0 0 1]) ;
  outstruct = pcaSM(mdatl,paramstruct) ;
  mpc = getfield(outstruct,'mpc') ;

  %  Set up Brushed Color Matrix
  %
  mcolor = ones(n,1) * [0 0 0] ;
      %  Start with all black and update
  if ig == 1 ;
    vflag = (mpc(2,:) < -10)' ; 
    mcolor(vflag,:) = ones(sum(vflag),1) * [1 0 0] ;
    vflag = (mpc(2,:) > 11)' ; 
    mcolor(vflag,:) = ones(sum(vflag),1) * [0 0 1] ;
  elseif ig == 2 ;
    vflag = (mpc(1,:) > 0)' ; 
    mcolor(vflag,:) = ones(sum(vflag),1) * [1 0 0] ;
    vflag = (mpc(2,:) > 4)' ; 
    mcolor(vflag,:) = ones(sum(vflag),1) * [0 0 1] ;
  elseif ig == 3 ;
    vflag = (mpc(1,:) > 70)' ; 
    mcolor(vflag,:) = ones(sum(vflag),1) * [1 0 0] ;
    vflag = (mpc(3,:) < -11)' ; 
    mcolor(vflag,:) = ones(sum(vflag),1) * [0 0 1] ;
  elseif ig == 4 ;
    vflag = (mpc(2,:) > 4)' ; 
    mcolor(vflag,:) = ones(sum(vflag),1) * [0 0 1] ;
    vflag = (mpc(1,:) > 60)' ; 
    mcolor(vflag,:) = ones(sum(vflag),1) * [1 0 0] ;
  end ;


  if ipart == 1 ;    %  Raw Data Curves, original scale


    figure(1) ;
    clf ;


    nbp = size(mdat,1) ;
    vibp = (1:nbp)' ;
    vaxh = axisSM(vibp) ;
    vaxv = axisSM(mdat) ;
    vaxv(1) = 0 ;

    plot(mdat,'-') ;
    title('Gene = CDKN2A') ;
    xlabel('exonic nt number, not genomic position') ;
    ylabel('RNA read depth') ;
    axis([vaxh(1) vaxh(2) vaxv(1) vaxv(2)]) ;

    orient landscape ;
    savestr = ['LungCancer2011ip' num2str(ipart) 'RawDatOverlay'] ;
    print('-dpsc2',savestr) ;


  elseif ipart == 2 ;    %  Raw Data Curves, log10 scale

    nbp = size(mdatl,1) ;
    vibp = (1:nbp)' ;
    vaxh = axisSM(vibp) ;
    vaxv = axisSM(mdatl) ;
    vaxv(1) = 0 ;
    plot(mdatl,'-') ;
    title('Gene = CDKN2A') ;
    xlabel('exonic nt number, not genomic position') ;
    ylabel('log_{10}(RNA read depth + 1)') ;
    axis([vaxh(1) vaxh(2) vaxv(1) vaxv(2)]) ;

    orient landscape ;
    savestr = ['LungCancer2011ip' num2str(ipart) 'Log10DatOverlay'] ;
    print('-dpsc2',savestr) ;


  elseif ipart == 3 ;    %  PCA scatterplot

    nbp = size(mdatl,1) ;
    titlecellstr = {{'Lung Cancer - Next Gen' ...
                     ['Gene = ' GeneNamesS{ig}] ...
                     ['n = ' num2str(n) 'patients'] ... 
                     ['d = ' num2str(nbp) ' base pairs']}} ;
    savestr = ['LungCancer2011ip' num2str(ipart) 'PCAScatPlot'] ;
    paramstruct = struct('npcadiradd',4, ...
                         'titlecellstr',titlecellstr, ...
                         'labelcellstr',{{'PC 1'; 'PC 2'; 'PC 3'; 'PC 4'}}, ...
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    scatplotSM(mdatl,[],paramstruct) ;


  elseif ipart == 4 ;    %  PCA scatterplot, brushed

    nbp = size(mdatl,1) ;
    titlecellstr = {{'Lung Cancer - Next Gen' ...
                     ['Gene = ' GeneNamesS{ig}] ...
                     ['n = ' num2str(n) 'patients'] ... 
                     ['d = ' num2str(nbp) ' base pairs']}} ;
    savestr = ['LungCancer2011ip' num2str(ipart) 'PCAScatPlotBrushed'] ;
    paramstruct = struct('npcadiradd',4, ...
                         'icolor',mcolor, ...
                         'titlecellstr',titlecellstr, ...
                         'labelcellstr',{{'PC 1'; 'PC 2'; 'PC 3'; 'PC 4'}}, ...
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    scatplotSM(mdatl,[],paramstruct) ;


  elseif ipart == 5 ;    %  Raw Data Curves, log10 scale, brushed

    %  Make data overlay - log10 Data
    %
    nbp = size(mdatl,1) ;
    vibp = (1:nbp)' ;
    vaxh = axisSM(vibp) ;
    vaxv = axisSM(mdatl) ;
    vaxv(1) = 0 ;

    hold on ;
    for i = 1:n ;
      if sum(mcolor(i,:)) > 0.5 ;
        plot(mdatl(:,i),'-','Color',mcolor(i,:)) ;
      else ;
        plot(mdatl(:,i),'k-') ;
      end ;
    end ;
    hold off ;
    title('Gene = CDKN2A') ;
    xlabel('exonic nt number, not genomic position') ;
    ylabel('log_{10}(RNA read depth + 1)') ;
    axis([vaxh(1) vaxh(2) vaxv(1) vaxv(2)]) ;

    orient landscape ;
    savestr = ['LungCancer2011ip' num2str(ipart) 'Log10DatOverlayBrushed'] ;
    print('-dpsc2',savestr) ;


  end ;    %  of inner ipart if-block


end ;   %  of outer ipart if-block





