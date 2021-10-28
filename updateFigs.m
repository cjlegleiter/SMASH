function updateFigs(matFile,scaleBarDist,scaleOffset,cBarOffset,legLabelOffset,...
    legBoxWidth)
%% Update figures based on existing outputs from the SMASH workflow
%
%% updateFigs.m:
% Helper function for producing final, updated figures based on output from the
% SMASH workflow. Input is a *.mat file from a previous SMASH run, along with
% some parameters to place a scale bar and legend to allow for tuning to
% different sites/dates.  Outputs are series of figures, each saved in an
% updatedFigs directory in *.fig, *.jpg, and *.eps format.
%
%% SYNTAX:
% updateFigs(matFile,scaleBarDist,scaleOffset,cBarOffset,legLabelOffset,legBoxWidth);
%
%% INPUTS:
% matFile:          String specifying the *.mat data file for the previous SMASH
%                   run for which updated figures will be produced
% scaleBarDist:     Scalar specifying the length of a scale bar to be placed on
%                   all maps, in units of m (so 10000 for 10 km)
% scaleOffset:      2 X 1 vector of x and y offsets for the scale bar relative
%                   to the default placement, will have to be set by trial and
%                   error for a particular image to avoid a poor placement
% cBarOffset:       Scalar with the offset of the colorbar from the classified
%                   map, will have to be set by trial and error for a particular
%                   image to avoid a poor placement (e.g., 0.05)
% legLabelOffset:   Scalar with the offset for legend labels on the classified
%                   map, will have to be set by trial and error for a particular
%                   image to avoid a poor placement (e.g., 0.02)
% legBoxWidth:      Scalar with the offset of the colorbar from the classified
%                   map, will have to be set by trial and error for a particular
%                   image to avoid a poor placement (e.g., 0.05)
%
%% OUTPUTS: a series of figures including
% RGB image
% Normalized difference chlorophyll index (NDCI) image
% Cyanobacterial index (CI) image
% Maximum RMSE sensitivity analysis
% MESMA-based classified map
% Histogram for distribution of genera from the MESMA-based classification
% Mask showing locations assigned to the genus selected by the user
% Fraction image for the genus selected by the user
% Fraction images for all genera and water
% False color composite based on two genera fractions and water fraction
% MESMA RMSE image
% Bar chart of RMSE by genus
%
%% NOTES:
%   This is a very high-level, post-processing helper function that takes output
%   from the SMASH workflow as input, so see preSMASH.m, runSMASH.m,
%   postSMASH.m, and tuneRmse.m for further detail.
%
%% FUNCTION SUMMARY:
% updateFigs(matFile,scaleBarDist,scaleOffset,cBarOffset,legLabelOffset,legBoxWidth);

%% CREDITS
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 10/20/2021
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\updateFigs.m


%% INTENT:
% Write a new function to read the existing figures, make a few changes, and
% save them as new .fig, .jpg, and .eps files. We can compile the .eps files in
% a series of LaTeX documents, one for each site/date, to include as supporting
% information for the manuscript. This approach will be cleaner and less labor
% intensive than the traditional PowerPoint figure dump.  Place the new figures
% into a separate "updateFigs" folder within each site/date folder we already
% have so that we can go back to the originals if necessary.

% Changes to make include:
    % Remove axes and add scale bar to all maps
    % Shift MESMA classification legend manually if necessary
    % Normalize class histogram by total number of pixels to get proportions
    % Replace the word "taxa" with "genera"
    % Italicize genera names
    % For taxa mask, modify title to "siteDateCode: Pixels classified as: Genera"
    % For taxa fractions, modify titke to "siteDateCode: Genera endmember fractions"
    % For all taxa fractions, modify overall title to "siteDateCode: Endmember fractions for each genus"
    % For taxaRGB, modify title to "siteDateCode: Endmember fraction composite ..."
    % For RMSE image, modify title to "siteDateCode: MESMA RMSE"
    % For mean RMSE by genus, make bars color-scaled based on proportion of
    % image clasified as that genus by making the color [0 1-proportion 0] for
    % each bar; also modify title to say genus rather than algal taxon
    % Also include tuneRMSE output plot.

%%    
%% Load the *.mat file with the data
% matFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810.mat';
load(matFile) %#ok<LOAD> 
siteDateCode=   [siteCode dateCode];

%%
%% Input (current figures) and output (updated figures) paths
newFigDir   =   [figDir(1:end-4) 'updateFigs']; %#ok<FNCOLND> 
mkdir(newFigDir)

%% 
%% Set scale bar length and shift relative to default position
% scaleBarDist=   10000;
% scaleOffset =   [5000 -5000];

%%
%% Update RGB image
% imgSubFile  =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810specSub.tif';
Hcube       =   hypercube(imgSubFile,wvlSub);
[rgb,iRgb]  =   colorize(Hcube,'Method','rgb','ContrastStretching',1);
rgbWvl      =   wvlSub(iRgb);
% Flip to account for georeferencing
rgb         =   flipud(rgb);
figure
imagesc(rgb,'XData',R.XWorldLimits,'YData',R.YWorldLimits,'alphadata',waterMask);
axis image; axis xy
title({[siteDateCode ' contrast-stretched RGB image'],...
       ['R: ' num2str(rgbWvl(1)) ' nm, G: ' num2str(rgbWvl(2)) ' nm, B: ' ...
       num2str(rgbWvl(3)) ' nm']});
% clear Hcube rgb iRgb

[hLine,hText]       =   scaleBar(scaleBarDist);
hLine.XData         =   hLine.XData + scaleOffset(1);
hLine.YData         =   hLine.YData + scaleOffset(2);
hText.Position(1)   =   hText.Position(1) + scaleOffset(1);
hText.Position(2)   =   hText.Position(2) + scaleOffset(2);
hText.String        =   [num2str(scaleBarDist/1000) ' km'];
axis off

set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%%
%% Update NDCI image
figure
imagesc(ndci,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
    'AlphaData',waterMask);
axis image; axis xy;
crameri('-bamako')
colorbar('location','eastoutside')
title([siteDateCode ': Normalized Difference Chlorophyll Index'])

[hLine,hText]       =   scaleBar(scaleBarDist);
hLine.XData         =   hLine.XData + scaleOffset(1);
hLine.YData         =   hLine.YData + scaleOffset(2);
hText.Position(1)   =   hText.Position(1) + scaleOffset(1);
hText.Position(2)   =   hText.Position(2) + scaleOffset(2);
hText.String        =   [num2str(scaleBarDist/1000) ' km'];
axis off

set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%%
%% Update CI image
figure
imagesc(ci,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
    'AlphaData',waterMask);
axis image; axis xy;
crameri('-bamako')
colorbar('location','eastoutside')
title([siteDateCode ': Cyanobacterial Index'])

[hLine,hText]       =   scaleBar(scaleBarDist);
hLine.XData         =   hLine.XData + scaleOffset(1);
hLine.YData         =   hLine.YData + scaleOffset(2);
hText.Position(1)   =   hText.Position(1) + scaleOffset(1);
hText.Position(2)   =   hText.Position(2) + scaleOffset(2);
hText.String        =   [num2str(scaleBarDist/1000) ' km'];
axis off

set(gcf,'name',[siteDateCode '_CI'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%%
%% Regnerate output figure from tuneRmse.m
load([matFile(1:end-4) 'rmseSensitivity.mat'],'maxRmse','classProp','Taxa')
if ~strcmp(Taxa{1},'Unclassified')
    Taxa    =   ['Unclassified'; Taxa];
end
TaxaLegend  =   Taxa;
for i = 2:length(Taxa)
    TaxaLegend{i}   =   ['{\it{' TaxaLegend{i} '}}'];
end

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
title([siteDateCode ': Maximum RMSE sensitivity analysis'])
legend(TaxaLegend)

set(gcf,'name',[siteDateCode '_maxRMSEsensitivity'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%%
%% Update classified map
figure
imagesc(modelSum,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & ~isnan(modelSum)); %#ok<*USENS> 
axis image; axis xy; axis off
title([siteDateCode ': MESMA-based classification']);   
cmap            =   crameri('bam',size(model,3));
cmap            =   [0 0 0; cmap];
colormap(cmap);
set(gca,'CLim',[0 size(model,3)]);

% Tweak placement of colorbar and legend
hCbar           =   colorbar('east');
hCbar.Ticks     =   [];
% cBarOffset      =   0.05;
hCbar.Position(1)=  hCbar.Position(1) + cBarOffset;
% A couple of hard-wired parameters for the legend layout
% legLabelOffset  =   0.02;
% legBoxWidth     =   0.05;
legBoxEdge      =   hCbar.Position(1)+legLabelOffset;
legBoxHeight    =   hCbar.Position(4)/(length(Taxa));
legBoxBot       =   hCbar.Position(2);
for i = 1:length(TaxaLegend)
    textBoxDim  =   [legBoxEdge legBoxBot+legBoxHeight*(i-1) ...
                     legBoxWidth legBoxHeight];
    h(i) = annotation('textbox',textBoxDim,'String',TaxaLegend(i),'FitBoxToText','on',...
               'fontname','arial','fontsize',12,'linestyle','none',...
               'verticalalignment','bottom','horizontalalignment','left');
end

[hLine,hText]       =   scaleBar(scaleBarDist);
hLine.XData         =   hLine.XData + scaleOffset(1);
hLine.YData         =   hLine.YData + scaleOffset(2);
hText.Position(1)   =   hText.Position(1) + scaleOffset(1);
hText.Position(2)   =   hText.Position(2) + scaleOffset(2);
hText.String        =   [num2str(scaleBarDist/1000) ' km'];
axis off

set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%%
%% Update class histogram
figure
taxaHist        =  histogram(modelSum(:),-0.5:length(Taxa)-0.5,...
                        'FaceColor',[0.4660 0.6740 0.1880]);
classProp       =   taxaHist.Values/sum(taxaHist.Values);
close
figure
taxaHist2       =  histogram('BinEdges',-0.5:length(Taxa)-0.5,'BinCounts',classProp,...
                        'FaceColor',[0.4660 0.6740 0.1880]); %#ok<NASGU> 
set(gca,'xlim',[-0.5 length(Taxa)-0.5])
set(gca,'xtick',0:length(Taxa)-1)
set(gca,'xticklabels',TaxaLegend)
ylabel('Proportion of classified image pixels')
title([siteDateCode ': MESMA-based distribution of genera']);

set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%%
%% Update selected taxa mask
figure
imagesc(model(:,:,taxaId),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & classMask);
crameri('-bamako');    
axis image; axis xy; axis off;
title([siteDateCode ': Pixels classified as ' TaxaLegend{taxaId+1}])

[hLine,hText]       =   scaleBar(scaleBarDist);
hLine.XData         =   hLine.XData + scaleOffset(1);
hLine.YData         =   hLine.YData + scaleOffset(2);
hText.Position(1)   =   hText.Position(1) + scaleOffset(1);
hText.Position(2)   =   hText.Position(2) + scaleOffset(2);
hText.String        =   [num2str(scaleBarDist/1000) ' km'];
axis off

set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%%
%% Update selected taxa fraction image
figure
imagesc(fractions(:,:,taxaId),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & ~isnan(fractions(:,:,taxaId)));
crameri('-bamako');    
axis image; axis xy; axis off;
set(gca,'clim',[0 1])
colorbar
title([siteDateCode ': ' TaxaLegend{taxaId+1} ' endmember fractions'])

[hLine,hText]       =   scaleBar(scaleBarDist);
hLine.XData         =   hLine.XData + scaleOffset(1);
hLine.YData         =   hLine.YData + scaleOffset(2);
hText.Position(1)   =   hText.Position(1) + scaleOffset(1);
hText.Position(2)   =   hText.Position(2) + scaleOffset(2);
hText.String        =   [num2str(scaleBarDist/1000) ' km'];
axis off

set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%%
%% Update fraction image mosaic for all genera
% For all taxa fractions, modify overall title to "siteDateCode: Endmember
% fractions for each genus"

% Display all fractions in a tiled layout, hardwired as 3 X 5 for now
figure
tiledlayout(3,5,'TileSpacing','none','Padding','compact');
titles  =   [TaxaLegend(2:end); {'Water'}];
% titles  =   [TaxaLegend; {'Water'}];
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
        hCbar.Position  =   [0.6 0.055 0.02 0.2853];
    end
end
% hTitle = title(t,[siteDateCode ': Endmember fractions for each genus'])
hTitle  =   text(1.5,0.9,...
                {[siteDateCode ':'],'Endmember fractions','for each genus'},...
                'Units','Normalized','Fontname','Arial','FontWeight','Bold',...
                'FontSize',13); %#ok<NASGU> 

set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%%
%% Update fraction false color composite
figure
imagesc(rgbFrac,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & (~isnan(rgbFrac(:,:,1)) | ~isnan(rgbFrac(:,:,2)) ...
        | ~isnan(rgbFrac(:,:,3))));
axis image; axis xy; axis off;
title([siteDateCode ': ' TaxaLegend{taxaId+1} ' endmember fractions'])
title({[siteDateCode ': Endmember fraction composite:'],...
    ['R: ' TaxaLegend{taxaIdR+1} ', G: ' TaxaLegend{taxaIdG+1} ', B: Water']})

[hLine,hText]       =   scaleBar(scaleBarDist);
hLine.XData         =   hLine.XData + scaleOffset(1);
hLine.YData         =   hLine.YData + scaleOffset(2);
hText.Position(1)   =   hText.Position(1) + scaleOffset(1);
hText.Position(2)   =   hText.Position(2) + scaleOffset(2);
hText.String        =   [num2str(scaleBarDist/1000) ' km'];
axis off

set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%% 
%% Update RMSE image
figure
imagesc(rmse,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & ~isnan(rmse));
axis image; axis xy; axis off; 
colormap(gray)
colorbar
title([siteDateCode ': MESMA RMSE'])

[hLine,hText]       =   scaleBar(scaleBarDist);
hLine.XData         =   hLine.XData + scaleOffset(1);
hLine.YData         =   hLine.YData + scaleOffset(2);
hText.Position(1)   =   hText.Position(1) + scaleOffset(1);
hText.Position(2)   =   hText.Position(2) + scaleOffset(2);
hText.String        =   [num2str(scaleBarDist/1000) ' km'];
axis off

set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%%
%% Update mean RMSE by genus
meanRmse    =   zeros(size(Taxa));
for i = 1:length(Taxa)
    inTaxa      =   find(modelSum == i-2 & waterMask);
    tmpRmse     =   rmse(inTaxa);
    meanRmse(i) =   mean(tmpRmse,'omitnan');
end
% For mean RMSE by genus, make bars color-scaled based on proportion of
% image clasified as that genus by making the color [0 1-proportion 0] for
% each bar; also modify title to say genus rather than algal taxon
figure
b   =   bar(meanRmse,'FaceColor','flat');
ylabel('RMSE')
set(gca,'xticklabels',TaxaLegend)
title([siteDateCode ': MESMA RMSE by genus'])
for i = 1:length(classProp)
    b.CData(i,:)    =   [0 classProp(i)/max(classProp) 0];
    if classProp(i) < 0.01
        b.CData(i,:)    =   [1 1 1];
    end
%     b.FaceAlpha(i)    =   classProp(i);   
end

set(gcf,'name',[siteDateCode '_RMSEbyTaxa'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

% This code resets the bar for any genus with a proportion less than 0.01 to
% white and makes bars for the more common genera more brightly green and the
% less common genera closer to black