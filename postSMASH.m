function [Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
          taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                      dateCode,waterMask,libraryFile,outputRoot)
% Process results from perform MESMA of algal taxa as part of HAB remote sensing workflow
%
%% postSMASH.m:
%   Import and process results from Multiple End Member Spectral Mixture
%   Analysis (MESMA) of algal taxa as part of a workflow for remote sensing of
%   Harmful Algal Blooms (HABs).  This function takes the outputs from
%   preSMASH.m and runSMASH.m as inputs and programmatically imports MESMA
%   results and creates summary figures.
%
%% SYNTAX:
%   [Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
%    taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
%                                      dateCode,waterMask,libraryFile,outputRoot);
%
%% INPUTS:
%   imgSubFile:     String with file name of output image used to perform MESMA
%   wvlSub:         Vector of wavelengths retained (units of nm)
%   R:              Geo-referencing object
%   siteCode:       String with two-letter site code
%   dateCode:       String with 8-digit date code (YYYYMMDD)
%   waterMask:      Binary water mask image
%   libraryFile:    String specifying spectral library file with algal taxa for
%                   use as end-members in MESMA; should *NOT* include water
%   outputRoot:     Root file name for MESMA outputs, with a time stamp used to
%                   identify the MESMA run
%
%% OUTPUTS:
%   Taxa:           Cell array of strings with names of algal taxa in the
%                   library that served as end members for MESMA
%   model:          MESMA output with the best model (number of bands = number
%                   of classes); each band contains the library spectra number
%                   per class
%   modelSum:       Single-band image summarizing the MESMA model image, with
%                   the value for each pixel representing to the class number;
%                   a classified image is also displayed
%   taxaHist:       Histogram object summarizing distribution (number of pixels)
%                   for each taxa in the MESMA-based classification; histogram
%                   is also displayed in a new figure
%   taxaId:         Index of the taxa for which classMask was produced, input by
%                   the user as an interactive prompt
%   classMask:      Binary mask highlighting pixels classified as the taxa
%                   specfied by the user via an interactive prompt
%   fractions:      Image with each band representing the MESMA fractions for
%                   each of the algal taxa, plus an additional (last) band for
%                   the water end member used as shade. A fraction image for the
%                   user-selected taxa is created along with a second figure
%                   showing fraction images for all end-members, including water
%   taxaIdR:        Index of the taxa displated as red in a color composite 
%                   fraction image, input by the user as an interactive prompt
%   taxaIdG:        Index of the taxa displated as green in a color composite 
%                   fraction image, input by the user as an interactive prompt
%   rgbFrac:        3D array of fractions for the two taxa selected for display
%                   as red and green in the first two bands and the water
%                   fraction as blue in the third band, displayed in a figure
%   rmse:           MESMA model Root Mean Squared Error image, displayed in a
%                   new figure
%   meanRmse:       An nTaxa X 1 vector of mean RMSE values for the pixels
%                   classified as each taxa, also displayed as a bar graph
%
%% NOTES:
% > Requires the MATLAB Image Processing and Mapping Toolboxes
% > The Hyperspectral Imaging Toolbox add-in must also be installed
% > Code contains comments with inputs used for Lake Owasco prototype
% > Takes output from preSMASH.m and runSMASH.m as input
% > For information on the core algorithm, see the Python MESMA documentation at
%   https://mesma.readthedocs.io/en/latest/userguide/mesma_cli.html
%
%% FUNCTION SUMMARY:
%   [Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
%    taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
%                                      dateCode,waterMask,libraryFile,outputRoot);

%% CREDITS:
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 06/18/2021
% 08/12/2021 - Updated to output the mean RMSE for the pixels classified as each
%              taxa
% Also see prototyping code in: C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SMASHowasco.m
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\postSMASH.m


%% Use Hyperspectral Image Processing Toolbox to import hypercube and display RGB
siteDateCode=   [siteCode dateCode];
Hcube       =   hypercube(imgSubFile,wvlSub);
[rgb,iRgb]  =   colorize(Hcube,'Method','rgb','ContrastStretching',1);
rgbWvl      =   wvlSub(iRgb);
% Flip to account for georeferencing
rgb         =   flipud(rgb);
figure
imagesc(rgb,'XData',R.XWorldLimits,'YData',R.YWorldLimits,'alphadata',waterMask);
axis equal; axis xy
xlabel('Easting (m)')
ylabel('Northing (m)')
title({[siteDateCode ' contrast-stretched RGB image'],...
       ['R: ' num2str(rgbWvl(1)) ' nm, G: ' num2str(rgbWvl(2)) ' nm, B: ' ...
       num2str(rgbWvl(3)) ' nm']});
clear Hcube   


%% Get list of taxa from library and display numeric codes
Taxa    =   readtable([libraryFile(1:end-3) 'csv']); 
Taxa    =   Taxa.Class;
% Allow for unclassified pixels
Taxa    =   [{'Unclassified'}; Taxa];
TaxaNums=   cell(size(Taxa));
for i = 1:length(Taxa)
    TaxaNums{i} =   [num2str(i-1) ': ' Taxa{i}];
end
disp('Algal taxa numeric codes and taxa names:')
disp(TaxaNums)
hFig    =   uifigure('Position',[50 600 250 300],'Name','Algal Taxa Listing');
uilistbox(hFig,'Items',TaxaNums,'Position',[5 5 240 290]);


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
% % Previous version
% % Reset unmodeled pixels to NaN;
% model(model==-1)=   NaN;
% New code on 10/5/2021
% Add one so that pixel values match class numbers and unclassified pixels have
% a value of zero
model           =   model + 1;
% Apply mask
model           =   applyMask(model,waterMask);


%% Display classified map with a handmade legend
% New code on 10/5/2021
modelSum        =   zeros(size(model,[1 2]));
for i = 1:size(model,3)
    tmp             =   model(:,:,i) > 0;
    modelSum(tmp)   =   i;
end
% Mask areas outside the water body as NaN's
modelSum(~waterMask | isnan(modelSum)) =   NaN;
% % Previous version
% modelSum                =   sum(model,3);
% modelSum(~waterMask)    =   NaN;
% % Add one so that pixel values match class numbers - done above on the model
% % array itself
% % modelSum                =   modelSum+1;
% % Except for the unclassified pixels, which will have a value of -1*nEm+1, so
% % -11 for our current library. Let's reset the unclassified pixels to zero.
% % iUnclass                =   abs(modelSum - (-1*length(Taxa)+2)) < 0.01;
% iUnclass                =   modelSum == -1*(length(Taxa)-2);
% modelSum(iUnclass)      =   0;
figure
imagesc(modelSum,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & ~isnan(modelSum));
axis equal; axis xy
xlabel('Easting (m)')
ylabel('Northing (m)')
title([siteDateCode ': MESMA-based classification']);   
impixelinfo    
% colormap(lines(size(model,3)));
% Update to allow for unclassified pixels represented by 0
% crameri('bam',size(model,3));
cmap            =   crameri('bam',size(model,3));
cmap            =   [0 0 0; cmap];
colormap(cmap);
set(gca,'CLim',[0 size(model,3)]);

%% Tweak placement of colorbar and legend
hCbar           =   colorbar('east');
hCbar.Ticks     =   [];
hCbar.Position(1)=  hCbar.Position(1)-0.15;
% A couple of hard-wired parameters for the legend layout
legLabelOffset  =   0.02;
legBoxWidth     =   0.05;
legBoxEdge      =   hCbar.Position(1)+legLabelOffset;
legBoxHeight    =   hCbar.Position(4)/(length(Taxa));
% legBoxBot       =   hCbar.Position(2);
legBoxBot       =   hCbar.Position(2)+0.0;
% Legend          =   Taxa;
Legend          =   TaxaNums;
for i = 1:length(Legend)
    textBoxDim  =   [legBoxEdge legBoxBot+legBoxHeight*(i-1) ...
                     legBoxWidth legBoxHeight];
    h(i) = annotation('textbox',textBoxDim,'String',Legend(i),'FitBoxToText','on',...
               'fontname','arial','fontsize',12,'linestyle','none',...
               'verticalalignment','bottom','horizontalalignment','left'); %#ok<AGROW,NASGU> 
end


%% Histogram of class assignments
figure
taxaHist =  histogram(modelSum(:),-0.5:length(Taxa)-0.5,'FaceColor',[0.4660 0.6740 0.1880]);
set(gca,'xlim',[-0.5 length(Taxa)-0.5])
set(gca,'xtick',0:length(Taxa)-1)
set(gca,'xticklabels',Taxa)
ylabel('Number of classified image pixels')
title([siteDateCode ': MESMA-based distribution of taxa']);


%% Specify taxa number to display area mapped as a single class
taxaId          =   input('Specify taxa ID number to display locations classified as that taxa: ');
% iTaxa           =   find(model(:,:,taxaId)==taxaId-1);
iTaxa           =   find(model(:,:,taxaId)==taxaId);
classMask       =   false(size(model(:,:,taxaId)));
classMask(iTaxa)=   true;
figure
imagesc(model(:,:,taxaId),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & classMask);
crameri('-bamako');    
axis image; axis xy; axis on;
xlabel('Easting (m)')
ylabel('Northing (m)')
% title(['MESMA pixels for which best model was for taxa: ' Taxa(taxaId)])
title(['MESMA pixels for which best model was for taxa: ' Taxa(taxaId+1)])


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


%% Use same taxa number to display and display fraction image for that taxa
figure
imagesc(fractions(:,:,taxaId),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & ~isnan(fractions(:,:,taxaId)));
axis image; axis xy; axis on; 
crameri('-bamako');
colorbar
xlabel('Easting (m)')
ylabel('Northing (m)')
% title(['MESMA fractions for taxa: ' Taxa(taxaId)])
title(['MESMA fractions for taxa: ' Taxa(taxaId+1)])


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
title(t,'MESMA fractions for each algal taxa')


%% After looking at all end member fractions (including water), select two taxa for RGB display
taxaIdR         =   input('Specify taxa ID number for fraction to display as red: ');
taxaIdG         =   input('Specify taxa ID number for fraction to display as green: ');
% Assume water (shade) is the last end member in the fractions stack
rgbFrac         =   fractions(:,:,[taxaIdR taxaIdG size(fractions,3)]);
figure
imagesc(rgbFrac,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & (~isnan(rgbFrac(:,:,1)) | ~isnan(rgbFrac(:,:,2)) ...
        | ~isnan(rgbFrac(:,:,3))));
axis image; axis xy; axis on; 
xlabel('Easting (m)')
ylabel('Northing (m)')
% title({'Fraction composite:',['R:' Taxa{taxaIdR} ' G:' Taxa{taxaIdG} ...
%        ' B:Water']})
title({'Fraction composite:',['R:' Taxa{taxaIdR+1} ' G:' Taxa{taxaIdG+1} ...
       ' B:Water']})

impixelinfo   


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


%% Display RMSE image
figure
imagesc(rmse,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & ~isnan(rmse));
axis image; axis xy; axis on; 
colormap(gray)
colorbar
xlabel('Easting (m)')
ylabel('Northing (m)')
title('MESMA RMSE')


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
title([siteDateCode ': MESMA RMSE by algal taxon'])