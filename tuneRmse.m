function [Taxa,MesmaOut,classProp] = ...
    tuneRmse(maxRmse,runDelay,dataDir,siteCode,dateCode,libraryFile,...
             classNameField,imgSubFile,shadeFile,imgScaleFactor,libScaleFactor,...
             minEmFrac,maxEmFrac,minShadeEmFrac,maxShadeEmFrac,waterMask,R)
% Tune MESMA maximum RMSE threshold to avoid false positives while retaining true positives
%
%% tuneRmse.m:
%   High-level function for tuning the maximum RMSE constraint used in MESMA to
%   avoid false positives (i.e., assigning a pixel to an algal taxon when that
%   taxon is not actually present) that could lead to misleading
%   classifications.  Conversely, we want to retain true positives (i.e.,
%   successfully detecting an algal taxon that is present).  This function helps
%   to strike a balance between a very small max RMSE threshold that would avoid
%   false positives by not classifying anything and a max RMSE threshold that is
%   too large and leads to pixels that don't really have any algae (or at least
%   not any included in the library) being assigned to a taxa that isn't a good
%   match to the observed image spectrum.
%
%% SYNTAX:
% [Taxa,MesmaOut,classProp] = ...
%     tuneRmse(maxRmse,runDelay,dataDir,siteCode,dateCode,libraryFile,...
%              classNameField,imgSubFile,shadeFile,imgScaleFactor,libScaleFactor,...
%              minEmFrac,maxEmFrac,minShadeEmFrac,maxShadeEmFrac,waterMask,R);
%
%% INPUTS: 
%   maxRmse:        Vector of maximum RMSE constraint values to be used in a
%                   sequence of MESMA runs, one for each value of maxRMSE
%   runDelay:       Scalar number of seconds to pause MATLAB execution while
%                   waiting for the MESMA algorithm to finish running in a
%                   command window
%   dataDir:        String specifying directory where outputs from this function
%                   will be created
%   siteCode:       String with two-letter site code
%   dateCode:       String with 8-digit date code (YYYYMMDD)
%   libraryFile:    String specifying spectral library file with algal taxa for
%                   use as end-members in MESMA; should *NOT* include water
%   classNameField: String specifying name used to identify end-members in
%                   library, usually 'Class'
%   imgSubFile:     String with file name of image to use in MESMA, output from
%                   preSMASH.m
%   shadeFile:      String specifying spectral library file with a single water 
%                   spectrum for use as the shade end-member in MESMA
%   imgScaleFactor: Scalar with scale factor for image, typically 1 or 10000
%   libScaleFactor: Scalar with scale factor for library, typically 1
%   minEmFrac:      Minimum end member fraction constraint used in MESMA model
%   maxEmFrac:      Maximum end member fraction constraint used in MESMA model
%   minShadeEmFrac: Minimum shade end member fraction constraint used in MESMA model
%   maxShadeEmFrac: Maximum shade end member fraction constraint used in MESMA model
%   waterMask:      Binary water mask image
%   R:              Geo-referencing object
%
%% OUTPUTS:
%   Taxa:           Cell array of strings with names of algal taxa in the
%                   library that served as end members for MESMA
%   MesmaOut:       Structure array of MESMA output for each max RMSE value in
%                   the following fields:
%   .model:         MESMA output with the best model (number of bands = number
%                   of classes); each band contains the library spectra number
%                   per class
%   .modelSum:      Single-band image summarizing the MESMA model image, with
%                   the value for each pixel representing to the class number
%   .taxaHist:      Histogram object summarizing distribution (number of pixels)
%                   for each taxa in the MESMA-based classification; histogram
%                   is also displayed in a new figure
%   .fractions:     Image with each band representing the MESMA fractions for
%                   each of the algal taxa, plus an additional (last) band for
%                   the water end member used as shade. A figure is created to 
%                   show fraction images for all end-members, including water
%   .rmse:          MESMA model Root Mean Squared Error image, displayed in a
%                   new figure
%   .meanRmse:      An nTaxa X 1 vector of mean RMSE values for the pixels
%                   classified as each taxa, also displayed as a bar graph
%   classProp:      nMaxRmse X nEndMember matrix with the proportion of the
%                   image pixels assigned to each taxon, or left unclassified;
%                   also plotted in a figure as a function of maxRmse
%
%% NOTES:
% This is a high-level function that calls the main routine runSMASH.m and uses
% much of the core code from postSMASH.m, so see the documentation for those
% functions for further detail.
%
%% FUNCTION SUMMARY:
% [Taxa,MesmaOut,classProp] = ...
%     tuneRmse(maxRmse,runDelay,dataDir,siteCode,dateCode,libraryFile,...
%              classNameField,imgSubFile,shadeFile,imgScaleFactor,libScaleFactor,...
%              minEmFrac,maxEmFrac,minShadeEmFrac,maxShadeEmFrac,waterMask,R);

%% CREDITS:
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 10/05/2021
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\tuneRmse.m


%% Detroit Lake test case
% load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813\DL20200813.mat')
% % Initially we used a max RMSE threshold of 0.025, so let's reduce this to
% % something very small like 0.01 and re-run the MESMA
% maxRmse     =   linspace(0.001,0.005,5);
% runDelay    =   10;

%%
%% Get list of taxa from library
Taxa    =   readtable([libraryFile(1:end-3) 'csv']); 
Taxa    =   Taxa.Class;
% Allow for unclassified pixels
Taxa    =   [{'Unclassified'}; Taxa];
TaxaNums=   cell(size(Taxa));
for i = 1:length(Taxa)
    TaxaNums{i} =   [num2str(i-1) ': ' Taxa{i}];
end

%%
%% Set siteDateCode
siteDateCode    =   [siteCode dateCode];

%%
%% loop over maxRmse values
classProp   =   zeros(length(maxRmse),length(Taxa));
for iRmse = 1:length(maxRmse)

    %% Run SMASH
    % Keep everything else the same and re-run the MESMA
    outputRoot  =   runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                              classNameField,imgSubFile,shadeFile,...
                              imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                              minShadeEmFrac,maxShadeEmFrac,maxRmse(iRmse));
    pause(runDelay)
    
    %%
    %% Strip down postSMASH.m to just give us what we need
    % [Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
    %           taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
    %                                       dateCode,waterMask,libraryFile,outputRoot)
    
    %%
    %% *MESMA file contains the best model for each pixel as library spectra number per class
    % This file has the best model [nb of bands = nb of classes] - each band
    % contains the library spectra number per class
    % Value of unmodeled pixels in output: -1
    % Value of pixels with no data in output: -2
    modelFile       =   outputRoot;
    model           =   envi2matlab(modelFile);
    % Flip to account for geo-referencing
    model           =   flipud(model);
    % Reset nodata pixels to NaN
    model(model==-2)=   NaN;
    % Add one so that pixel values match class numbers and unclassified pixels have
    % a value of zero
    model           =   model + 1;
    % Apply mask
    model           =   applyMask(model,waterMask);
    
    
    %% Summarize model with a classification
    modelSum        =   zeros(size(model,[1 2]));
    for i = 1:size(model,3)
        tmp             =   model(:,:,i) > 0;
        modelSum(tmp)   =   i;
    end
    % Mask areas outside the water body as NaN's
    modelSum(~waterMask | isnan(modelSum)) =   NaN;
    
    %%
    %% Histogram of class assignments
    figure
    taxaHist =  histogram(modelSum(:),-0.5:length(Taxa)-0.5,'FaceColor',[0.4660 0.6740 0.1880]);
    set(gca,'xlim',[-0.5 length(Taxa)-0.5])
    set(gca,'xtick',0:length(Taxa)-1)
    set(gca,'xticklabels',Taxa)
    ylabel('Number of classified image pixels')
    title({[siteDateCode ': MESMA-based distribution of taxa'], ...
           ['Maximum RMSE constraint: ' num2str(maxRmse(iRmse))]});
    
    %%
    %% *MESMA_fractions file contains the model's fractions, including shade
    % The model’s fractions [nb of bands = nb of classes + 1], including a shade fraction
    % Value of unmodeled pixels in output: 0
    % Value of pixels with no data in output: 0
    % fractionsFile   =   fullfile(dataDir,[siteDateCode 'SpecSubMESMA_fractions']);
    fractionsFile   =   fullfile([outputRoot '_fractions']);
    fractions       =   envi2matlab(fractionsFile);
    % Flip to account for geo-referencing
    fractions       =   flipud(fractions);
    % Reset nodata pixels to NaN
    fractions(fractions==0)=   NaN;
    % Reset unmodeled pixels to NaN;
    fractions(fractions==0)=   NaN;
    % Apply mask
    fractions       =   applyMask(fractions,waterMask);
    
    %%
    %% Display all fractions in a tiled layout, hardwired as 3 X 5 for now
    figure
    t       =   tiledlayout(3,5,'TileSpacing','none','Padding','compact');
    % titles  =   [Taxa; {'Water'}];
    titles  =   [Taxa(2:end); {'Water'}];
    for i = 1:size(fractions,3)
        nexttile
        imagesc(fractions(:,:,i),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
            'alphadata',waterMask & ~isnan(fractions(:,:,i)));
        axis image; axis xy; axis off; 
        set(gca,'clim',[0 1])
        crameri('-bamako');
        title(titles(i))
        if i == size(fractions,3)
            hCbar           =   colorbar('EastOutside');
            hCbar.Position  =   [0.8 0.055 0.025 0.2853];
        end
    end
    title(t,{'MESMA fractions for each algal taxa',...
             ['Maximum RMSE constraint: ' num2str(maxRmse(iRmse))]});
    
    %%
    %% *MESMA_rmse file contains the model's RMSE
    % The model’s RMSE [nb of bands = 1]
    % Value of unmodeled pixels in output: 9999
    % Value of pixels with no data in output: 9998
    rmseFile    =   [outputRoot '_rmse'];
    rmse        =   envi2matlab(rmseFile);
    % Flip to account for geo-referencing
    rmse        =   flipud(rmse);
    % Reset nodata pixels to NaN
    rmse(rmse==9998)=   NaN;
    % Reset unmodeled pixels to NaN;
    rmse(rmse==9999)=   NaN;
    % Apply mask
    rmse        =   applyMask(rmse,waterMask);
    
    %%
    %% Display RMSE image
    figure
    imagesc(rmse,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
            'alphadata',waterMask & ~isnan(rmse));
    axis image; axis xy; axis on; 
    colormap(gray)
    colorbar
    xlabel('Easting (m)')
    ylabel('Northing (m)')
    title({'MESMA RMSE',...
          ['Maximum RMSE constraint: ' num2str(maxRmse(iRmse))]});
    
    %%
    %% Calculate mean RMSE for pixels classified as each taxa
    meanRmse    =   zeros(size(Taxa));
    for i = 1:length(Taxa)
    %    inTaxa       =   find(modelSum == i & waterMask);
        inTaxa      =   find(modelSum == i-2 & waterMask);
        tmpRmse     =   rmse(inTaxa); %#ok<*FNDSB>
        meanRmse(i) =   mean(tmpRmse,'omitnan');
    end
    figure
    bar(meanRmse);
    ylabel('RMSE')
    set(gca,'xticklabels',Taxa)
    title({[siteDateCode ': MESMA RMSE by algal taxon'],...
           ['Maximum RMSE constraint: ' num2str(maxRmse(iRmse))]});
    
    %%
    %% Calculate proportion of pixels assigned to each taxa, or left unclassified
    classProp(iRmse,:)  =   taxaHist.Values/sum(taxaHist.Values);

    %% Allocate MESMA results to structure array
    MesmaOut(iRmse).model    =  model; %#ok<*AGROW> 
    MesmaOut(iRmse).modelSum =  modelSum;
    MesmaOut(iRmse).taxaHist =  taxaHist;
    MesmaOut(iRmse).fractions=  fractions;
    MesmaOut(iRmse).rmse     =  rmse;
    MesmaOut(iRmse).meanRmse =  meanRmse;

end


%%
%% Plot class proportions vs. maxRMSE threshold
figure
h       =   plot(maxRmse,classProp,'linewidth',2);
h(1).LineStyle  = '-';
h(2).LineStyle  = '--';
h(3).LineStyle  = ':';
h(4).LineStyle  = '-.';
h(5).LineStyle  = '-';
h(6).LineStyle  = '--';
h(7).LineStyle  = ':';
h(8).LineStyle  = '-.';
h(9).LineStyle  = '-';
h(10).LineStyle = '--';
h(11).LineStyle = ':';
h(12).LineStyle = '-.';
h(13).LineStyle = '-';
xlabel('Maximum RMSE threshold')
ylabel('Classified taxa proportions')
title([siteDateCode ': Sensitivity to maximum RMSE theshold'])
legend(TaxaNums)
