%% SMASH workflow for Owasco Lake example
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 05/26/2021
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SMASHowasco.m


%% Set data directory, two-letter site code, date, EPSG projection, and bands to keep
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco';
siteCode    =   'OW';
dateCode    =   '20200823';
% WGS 84 UTM Zone 18N is EPSG 32618;
epsgCode    =   32618;
siteDateCode=   [siteCode dateCode];
% Specify bands to keep to match convolved library
iKeepWvl    =   1:158;


%% Load original DESIS Geotiff image and crop to spatial subset - Faster in ENVI
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\rawImage\DESIS-HSI-L2A-DT0489826624_003-20200823T165204-V0210-SPECTRAL_IMAGE.tif';
% [img,R]     =   readgeoraster(imgFile);
% Loading the image this way is very slow,maybe try multibandread instead?
disp('Open this file in ENVI and use the ROI tool to establish spatial subset:')
disp(imgFile)
disp(['Save the ROI as: ' fullfile(dataDir,[siteDateCode '.roi'])])
% Open file in ENVI and then create an ROI for the subset of interest:
%   C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OwascoSubset.roi
% To do so, load the image to a display in ENVI, then from the display menu,
% select Tools-->Region Of Interest-->ROI tool. In the ROI tool, select the
% image as the window to create the roi and change the roi type to rectangle.

disp('Use ENVI to create a new spatial subset GEOTIFF based on this ROI.')
disp(['Save the spatial subset image as: ' fullfile(dataDir,[siteDateCode 'sub.tif'])])
% Used this ROI to export the data as a new GeoTiff. From the ENVI file menu,
% File-->Save File As-->TIFF/GeoTIFF, and in the resulting dialog select the
% original image, then click spatial subset, and point to the ROI you created.


%% Load the subset image file created in ENVI
imgSubFile  =   fullfile(dataDir,[siteDateCode 'sub.tif']);
[img,R]     =   readgeoraster(imgSubFile);
% Flip the image to make the geo-referencing correct
img         =   flipud(img);


%% Load wavelength information from ENVI header for the original file
% doc readEnviHeader
% Don't keep all of the output arguments since we'll be using the subset image,
% we really just need the band information
% [samples,lines,bands,offset,dataType,interleave,byteOrder,MapInfo,wvl] = ...
%         readEnviHeader(imgFile);
[~,~,~,~,~,~,~,~,wvl]   =   readEnviHeader(imgFile);


%% Calculate NDWI and display NDWI image
% (Band 63 - band 182) / (band 63 + band 182), where 63 is green (560.6 nm) and
% 182 is NIR (865.4 nm)
iGrnWvl =   63;
iNirWvl =   182;
img     =   double(img);
ndwi    =   (img(:,:,iGrnWvl) - img(:,:,iNirWvl))./(img(:,:,iGrnWvl) + img(:,:,iNirWvl));
figure
imshow(ndwi,'XData',R.XWorldLimits,'YData',R.YWorldLimits,'DisplayRange',[]);
axis xy; axis on; colorbar
title([siteDateCode ': Normalized Difference Water Index'])
xlabel('Easting (m)')
ylabel('Northing (m)')


%% Create water mask based on NDWI and apply to image
% doc buildWaterMask1band
% [waterMask,bandThresh] = buildWaterMask1band(img,imXdata,imYdata,wvl)
% Note that we have to multiply by 1000 to get a reasonable data range
[waterMask,bandThresh] = buildWaterMask1band(ndwi*1000,R.XWorldLimits,R.YWorldLimits,1);
% Used [-500 1000]
% Apply mask to image
imgMask =   applyMask(img,waterMask);


%% Display the masked image for a green band
figure
imagesc(imgMask(:,:,iGrnWvl),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
    'AlphaData',waterMask);
axis image; axis xy; axis on;
cmap    =   crameri('-bamako');
crameri('-bamako')
colorbar
title([siteDateCode ': Masked image for ' num2str(wvl(iGrnWvl)) ' nm band'])
xlabel('Easting (m)')
ylabel('Northing (m)')


%% Calculate NDCI and display result
% (Band 120 - Band 104) / (Band 120 + Band 104)
% Where 120 is NIR (706.7 nm) and 104 is Red (665.3 nm)
iRedWvl =   104;
iNirWvl2=   120;
ndci    =   (imgMask(:,:,iNirWvl2) - imgMask(:,:,iRedWvl))./...
            (imgMask(:,:,iNirWvl2) + imgMask(:,:,iRedWvl));
ndci    =   applyMask(ndci,waterMask);
figure
imagesc(ndci,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
    'AlphaData',waterMask);
axis image; axis xy; axis on;
cmap    =   crameri('-bamako');
crameri('-bamako')
colorbar
title([siteDateCode ': Normalized Difference Chlorophyll Index'])
xlabel('Easting (m)')
ylabel('Northing (m)')


%% Export mask, NDWI, and NDCI as Geo-TIFF's
% Note that we have to reflip so that the exported images are in the correct
% orientation using the same R geo-referencing structure created when we
% imported the image in the first place
% Export mask as a geotiff
maskFile=   fullfile(dataDir,[siteDateCode 'mask.tif']);
geotiffwrite(maskFile,flipud(waterMask),R,'CoordRefSysCode',epsgCode);
% Export NDWI image
ndwiFile=   fullfile(dataDir,[siteDateCode 'ndwi.tif']);
geotiffwrite(ndwiFile,flipud(ndwi),R,'CoordRefSysCode',epsgCode);
% Export NDCI
ndciFile=   fullfile(dataDir,[siteDateCode 'ndci.tif']);
geotiffwrite(ndciFile,flipud(ndci),R,'CoordRefSysCode',epsgCode);


%% Spectral subset to match convolved library and export
imgSpecSub  =   img(:,:,iKeepWvl);
imgSpecSub  =   applyMask(imgSpecSub,waterMask);
wvlSub      =   wvl(iKeepWvl);
imgSubFile  =   fullfile(dataDir,[siteDateCode 'SpecSub.tif']);
% Note that we have to reflip so that the exported images are in the correct
% orientation using the same R geo-referencing structure created when we
% imported the image in the first place
geotiffwrite(imgSubFile,flipud(imgSpecSub),R,'CoordRefSysCode',epsgCode);
% Export a simple text file with wavelengths we can import to ENVI
writematrix(wvlSub,fullfile(dataDir,'wvlSubList.txt'));


%% Now on to MESMA in QGIS
% https://mesma.readthedocs.io/en/latest/
% Note that the MESMA tool can be run using a settings file, which we might be
% able to generate programmatically to minimize GUI interaction (and the use of
% QGIS at all). For an example of this type of file, see:
%   C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\MesmaInput2.txt


%% Set parameters we will pass to the MESMA settings text file read by QGIS
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\AlgaeTaxaOnly.sli';
% Metadata field with end-member names
classNameField  =   'Class';
% Constraints
%   Minimum end member fraction
minEmFrac       =   -0.05;
%   Maximum end member fraction
maxEmFrac       =   1.05;
%   Minimum shade end member fraction
minShadeEmFrac  =   -0.05;
%   Maximum shade end member fraction
maxShadeEmFrac  =   1.05;
%   Maximum RMSE
maxRMSE         =   0.025;
% Name of shade end member file
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\DetroitLakeWaterOnly.sli';


%% Information on QGIS MESMA tool API - copied here for convenience:

% See the website: https://mesma.readthedocs.io/en/latest/api/mesma.html

% And source Python code:
% C:\Users\cjl\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\mesma\core\mesma.py

% Execute MESMA. Process input and output.
% 
% In case band weighing or band selection algorithms are used, no residual image
% or residual constraints can be used.
% 
% Returns 3 images [a pixels * b pixels * c bands]:
% 
% the best model [nb of bands = nb of classes] - each band contains the library spectra number per class
% the model’s fractions [nb of bands = nb of classes + 1], including a shade fraction
% the model’s RMSE [nb of bands = 1]
% [optional] a residual image

% Value of unmodeled pixels in output:
% models: -1
% fractions: 0
% rmse: 9999
% residual_image: 0

% Value of pixels with no data in output:
% models: -2
% fractions: 0
% rmse: 9998
% residual_image: 0

% Parameters:	
% image – image, scaled to reflectance, without bad bands
% library – spectral library with spectra as columns, scaled to reflectance, without bad bands
% look_up_table – all endmember combinations (=models) for MESMA; ordered per complexity level and per class-model; n_models x n_endmembers
% em_per_class – a list of all library indices per endmember class
% constraints – min + max endmember fraction, min + max shade fraction, max rmse, residual reflectance threshold + max number of consecutive bands exceeding threshold. set value to -9999 if not used.
% no_data_pixels – indices of pixels that contain no data (result of np.where)
% shade_spectrum – single spectrum of photometric shade
% fusion_value – only select a model of higher complexity (e.g. 3-EM over 2-EM) of the RMSE is better with at least this value
% residual_image – output the residuals as an image (ignored when using band weighing or -selection)
% use_band_weighing – use the weighted linear spectral mixture analysis (Somers et al, 2009)
% use_band_selection – use the bands selection algorithm (Somers et al, 2010)
% bands_selection_values – correlation threshold and decrease for the band selection algorithm
% log – log function
% Returns:	
% images with the best model for each pixel, the model fractions and rmse belonging {+ evt. residuals)


%% Template for programmatic generation of MESMA setttings file for use in QGIS
% See the example: C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\MesmaInput2.txt
MesmaSettings   =   char(...
'MESMA Settings: ',...
'   ',... 
libraryFile,...% Name of input library
classNameField,...% Metadata field with end-member names
[num2str(minEmFrac) ', ' num2str(maxEmFrac) ', ' num2str(minShadeEmFrac) ', ' ...
num2str(maxShadeEmFrac) ', ' num2str(maxRMSE) ', -9999, -9999'],...% List of constraints:  min + max endmember fraction, min + max shade fraction, max rmse, residual reflectance threshold + max number of consecutive bands exceeding threshold. set value to -9999 if not used.
fullfile(dataDir,[siteDateCode 'SpecSub.tif']),...% Name of input image
fullfile(dataDir,[siteDateCode 'SpecSubMESMA']),...% Root name of MESMA output files
'True',...% Logical flag to open result in QGIS
'1',...% Spectral library reflectance scale factor
'0.007',...% Fusion threshold (not used)
shadeFile,...% Try name of non-photometric shade here
'1',...% Reflectance scale factor for shade end-member
'10000',...% Image reflectance scale factor
'False',...% Logical flag to create residual image
'False',...% Logical flag to use spectral weighting
'False',...% Logical flag to use band selection (Stable Zone Unmixing)
'0.99',...% Correlation threshold for Stable Zone Unmixing
'0.01',...% Correlation decrease for Stable Zone Unmixing
'x-EM Levels:',...% The rest sets up which models and end-members will be used and shouldn't change as long as the number of end-members doesn't change and we don't decide to exclude some end members on a case by case basis
'False, False, True, True, False, False, False, False, False, False, False, False, False, False, False',...
'classes per x-EM Level:',...
'2-EMs: True, True, True, True, True, True, True, True, True, True, True, True, True',...
'3-EMs: True, True, True, True, True, True, True, True, True, True, True, True, True',...
'4-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False',...
'5-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False',...
'6-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False',...
'7-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False',...
'8-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False',...
'9-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False',...
'10-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False',...
'11-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False',...
'12-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False',...
'13-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False',...
'14-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False',...
'models per x-EM Level:',...
'2-EMs: True, True, True, True, True, True, True, True, True, True, True, True, True',...
'3-EMs: True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True',...
'4-EMs: ',...
'5-EMs: ',...
'6-EMs: ',...
'7-EMs: ',...
'8-EMs: ',...
'9-EMs: ',...
'10-EMs: ',...
'11-EMs: ',...
'12-EMs: ',...
'13-EMs: ',...
'14-EMs:     ');
MesmaSettings   =   cellstr(MesmaSettings);
% writecell(MesmaSettings,'tmp.txt')
fid             =   fopen(fullfile(dataDir,[siteDateCode 'MesmaSettings.txt']),'w');
formatSpec      =   '%s\n';
for i = 1:length(MesmaSettings)
    fprintf(fid,formatSpec,MesmaSettings{i});
end
fclose(fid);
% type tmp.txt
% Looks like this worked and can be used to do a run in QGIS








%% **** 6/1/2021: Pick up with import of fraction images and other outputs from MESMA


%% Use Hyperspectral Image Processing Toolbox to import hypercube and display RGB
Hcube       =   hypercube(fullfile(dataDir,[siteDateCode 'SpecSub.tif']),wvlSub);
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


%% Get list of taxa from library and display numeric codes
Taxa    =   readtable([libraryFile(1:end-3) 'csv']);
Taxa    =   Taxa.Class;
TaxaNums=   cell(size(Taxa));
for i = 1:length(Taxa)
    TaxaNums{i} =   [num2str(i) ': ' Taxa{i}];
end
disp('Algal taxa numeric codes and taxa names:')
disp(TaxaNums)
hFig    =   uifigure('Position',[50 600 250 250],'Name','Algal Taxa Listing');
hLabel  =   uilistbox(hFig,'Items',TaxaNums,'Position',[5 5 240 240]);


%% *MESMA file contains the best model for each pixel as library spectra number per class
% This file has the best model [nb of bands = nb of classes] - each band
% contains the library spectra number per class
% Value of unmodeled pixels in output: -1
% Value of pixels with no data in output: -2
modelFile       =   fullfile(dataDir,[siteDateCode 'SpecSubMESMA']);
model           =   envi2matlab(modelFile);
% Flip to account for geo-referencing
model           =   flipud(model);
% Reset nodata pixels to NaN
model(model==-2)=   NaN;
% % Reset unmodeled pixels to NaN;
% model(model==-1)=   NaN;
% Apply mask
model           =   applyMask(model,waterMask);


%% Display classified map with a handmade legend
modelSum    =   zeros(size(model,[1 2]));
for i = 1:size(model,3)
    iTaxa       =   find(model(:,:,i)==i-1);
    tmp         =   zeros(size(model,[1 2]));
    tmp(iTaxa)  =   i; %#ok<FNDSB>
    modelSum    =   modelSum + tmp;
end
modelSum(~waterMask)    =   NaN;
figure
imagesc(modelSum,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask | modelSum==0);
map             =   colormap(lines(size(model,3)+1));
hCbar           =   colorbar;
hCbar.Ticks     =   [];
% A couple of hard-wired parameters for the legend layout
legLabelOffset  =   -0.1;
legBoxWidth     =   0.05;
legBoxEdge      =   hCbar.Position(1)+legLabelOffset;
legBoxHeight    =   hCbar.Position(4)/(length(Taxa)+1);
legBoxBot       =   hCbar.Position(2);
Legend          =   [{'Unclassified'}; Taxa];
for i = 1:length(Legend)
    textBoxDim  =   [legBoxEdge legBoxBot+legBoxHeight*(i-1) ...
                     legBoxWidth legBoxHeight];
    annotation('textbox',textBoxDim,'String',Legend(i),'FitBoxToText','on',...
               'fontname','arial','fontsize',12,'linestyle','none',...
               'verticalalignment','bottom','horizontalalignment','left');
end
axis equal; axis xy
xlabel('Easting (m)')
ylabel('Northing (m)')
title([siteDateCode ': MESMA-based classification']);


%% Specify taxa number to display area mapped as a single class
taxaId          =   input('Specify taxa ID number to display: ');
iTaxa           =   find(model(:,:,taxaId)==taxaId-1);
classMask       =   false(size(model(:,:,taxaId)));
classMask(iTaxa)=   true;
figure
imagesc(model(:,:,taxaId),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & classMask);
crameri('-bamako')    
axis image; axis xy; axis on;
xlabel('Easting (m)')
ylabel('Northing (m)')
title(['MESMA pixels with best model for taxa: ' Taxa(taxaId)])


%% *MESMA_fractions file contains the model's fractions, including shade
% The model’s fractions [nb of bands = nb of classes + 1], including a shade fraction
% Value of unmodeled pixels in output: 0
% Value of pixels with no data in output: 0
fractionsFile   =   fullfile(dataDir,[siteDateCode 'SpecSubMESMA_fractions']);
fractions       =   envi2matlab(fractionsFile);
% Flip to account for geo-referencing
fractions       =   flipud(fractions);
% Reset nodata pixels to NaN
fractions(fractions==0)=   NaN;
% % Reset unmodeled pixels to NaN;
% fractions(fractions==0)=   NaN;
% Apply mask
fractions       =   applyMask(fractions,waterMask);


%% Specify taxa number to display and display fraction image for that taxa
taxaId          =   9;
figure
imagesc(fractions(:,:,taxaId),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask);
axis image; axis xy; axis on; 
crameri('-bamako');
colorbar
xlabel('Easting (m)')
ylabel('Northing (m)')
title(['MESMA fractions for taxa: ' Taxa(taxaId)])


%% *MESMA_rmse file contains the model's RMSE
% The model’s RMSE [nb of bands = 1]
% Value of unmodeled pixels in output: 9999
% Value of pixels with no data in output: 9998
rmseFile   =   fullfile(dataDir,[siteDateCode 'SpecSubMESMA_rmse']);
rmse       =   envi2matlab(rmseFile);
% Flip to account for geo-referencing
rmse       =   flipud(rmse);
% Reset nodata pixels to NaN
rmse(rmse==9998)=   NaN;
% Reset unmodeled pixels to NaN;
rmse(rmse==9999)=   NaN;
% Apply mask
rmse       =   applyMask(rmse,waterMask);


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


%% Save what we have so far
clear c cmap ans fid i img lat lon p proj r wvl
save SMASHowasco.mat


%% ***** 6/3/2021: Try using only models that include water 
% Use the updated settings file template ForcedWater.txt, then experiment with
% constraints on shade fraction and RMSE


%% Simulate mixtures to see if we can reproduce known fractions


%% Add water to algal spectral library and save as an image file in ENVI, and import
% First we need to combine the algal taxa and the water end member from Detroit
% Lake into a single spectral library in ENVI and then save this new merged
% spectral library as an image file that we can read into MATLAB
% doc hypercube
mergedLibraryFile   =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\Algae+WaterLibraryAsImage.dat';
library             =   envi2matlab(mergedLibraryFile);
% This is a 14 row by 158 column array with each end member as a row and each
% wavelength as a column. The order of the taxa should be the same as before,
% but with water as the 14th end member now:
%     {'1: Anabaena'                  }
%     {'2: Aphanizomenon'             }
%     {'3: Cylindrospermum-corrrected'}
%     {'4: Eucapsis'                  }
%     {'5: Gleocapsa'                 }
%     {'6: Gleotrichia'               }
%     {'7: Lyngbya'                   }
%     {'8: Merismopedia'              }
%     {'9: Microcystis'               }
%     {'10: Nostoc'                   }
%     {'11: Oscillatoria'             }
%     {'12: Spirulina'                }
%     {'13: Tolypothrix'              }
% 14 is water
figure
plot(wvlSub,library)
xlabel('Wavelength (nm)')
ylabel('Reflectance (%)')
title('Merged algae and water spectral library')
legend([TaxaNums; {'14: Detroit Lake water'}])


%% Generate mixtures with known fractions
iTaxa2mix1  =   9;  % Microcystis
iTaxa2mix2  =   5;  % Gleocapsa
iWater      =   14;
fracTaxa1   =   0.5;
fracTaxa2   =   0.4;
fracWater   =   1-fracTaxa1-fracTaxa2;
mix         =   library(iTaxa2mix1,:)*fracTaxa1 + ...
                library(iTaxa2mix2,:)*fracTaxa2 + ...
                library(iWater,:)*fracWater;
hold on
plot(wvlSub,mix,'linewidth',2)
% Looks like this is working OK


%% Set up an image with known simulated mixtures
% Use the same dimensions as our existing image so we can keep the
% geo-referencing info etc.: 580 X 285
% 50 rows each for a water fraction of 0 to 1 in steps of 0.1
fracWater   =   (0:0.1:1)';
% Fractions for the two algae evenly spaced between the two end members for the
% non-water portion of each pixel
fracTaxa1   =   linspace(0,1,280);
fracTaxa2   =   1-fracTaxa1;
% Now allocate the water fractions to an array with the same dimensions as our
% actual image
fracWaterImg=   zeros(size(imgSpecSub,[1 2]));
for i = 1:length(fracWater)
    fracWaterImg(i+50*(i-1):i+50*i,:)   =   fracWater(i);
end
% Clean up edge of water fraction image to have strip of zeros like the taxa
% fraction images will
fracWaterImg(:,280:end) =   0;
figure
subplot(1,3,1)
imshow(fracWaterImg)
colorbar
title('Water Fraction')
% Same kind of thing for the taxa fractions, adjusted to account for the water
% fraction already allocated for that row
fracTaxa1img=   zeros(size(imgSpecSub,[1 2]));
fracTaxa2img=   zeros(size(imgSpecSub,[1 2]));
for i = 1:length(fracWater)
    fracTaxa1img(i+50*(i-1):i+50*i,1:length(fracTaxa1)) =   ...
        repmat(fracTaxa1,51,1)*(1-fracWater(i));
    fracTaxa2img(i+50*(i-1):i+50*i,1:length(fracTaxa2)) =   ...
        repmat(fracTaxa2,51,1)*(1-fracWater(i));    
end
subplot(1,3,2)
imshow(fracTaxa1img)
title('Taxa #1 Fraction')
colorbar
subplot(1,3,3)
imshow(fracTaxa2img)
title('Taxa #2 Fraction')
colorbar
figure
imshow(fracWaterImg+fracTaxa1img+fracTaxa2img)
colorbar
% Looks good, so next we make the mixtures


%% Simulate mixtures based on the predefined water and two taxa fractions
% Set up loop over an array of the same dimensions as our image
simSpec     =   zeros(size(imgSpecSub));
for i = 1:size(imgSpecSub,1)
   for j = 1:size(imgSpecSub,2)
       mixPix           =   library(iTaxa2mix1,:)*fracTaxa1img(i,j) + ...
                            library(iTaxa2mix2,:)*fracTaxa2img(i,j) + ...
                            library(iWater,:)*fracWaterImg(i,j);
       simSpec(i,j,:)   =   mixPix;
   end
end
% figure
% imagesc(flipud(simSpec(:,:,69)),'XData',R.XWorldLimits,'YData',R.YWorldLimits);
% axis image; axis xy
% colorbar
% Looks OK
simCube     =   hypercube(flipud(simSpec),wvlSub);
simRgb      =   colorize(simCube,"Method","rgb","ContrastStretching",false);
imshow(simRgb,'XData',R.XWorldLimits,'YData',R.YWorldLimits);
axis image; axis xy


%% Next write this simulated image out as a GEOTIFF we can use for MESMA
% Note that we might have to flip so that the exported image is in the correct
% orientation using the same R geo-referencing structure created when we
% imported the actual image in the first place
simFile     =   fullfile(dataDir,[siteDateCode 'simSpec.tif']);
geotiffwrite(simFile,simSpec,R,'CoordRefSysCode',epsgCode);
% Also export the fraction images for the two taxa and water
fracTax1file=   fullfile(dataDir,[siteDateCode 'fracTaxa1.tif']);
geotiffwrite(fracTax1file,fracTaxa1img,R,'CoordRefSysCode',epsgCode);
fracTax2file=   fullfile(dataDir,[siteDateCode 'fracTaxa2.tif']);
geotiffwrite(fracTax2file,fracTaxa2img,R,'CoordRefSysCode',epsgCode);
fracWatrFile=   fullfile(dataDir,[siteDateCode 'fracWater.tif']);
geotiffwrite(fracWatrFile,fracWaterImg,R,'CoordRefSysCode',epsgCode);



%% Now set up MESMA run for the simulated image, with forced water models


%% Set parameters we will pass to the MESMA settings text file read by QGIS
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\Algae+Water.sli';
% Metadata field with end-member names
classNameField  =   'Class';
% Constraints
%   Minimum end member fraction
minEmFrac       =   -0.05;
%   Maximum end member fraction
maxEmFrac       =   1.05;
%   Minimum shade end member fraction
minShadeEmFrac  =   0;
%   Maximum shade end member fraction
maxShadeEmFrac  =   0.01;
%   Maximum RMSE
maxRMSE         =   0.025;
% Name of shade end member file
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\DetroitLakeWaterOnly.sli';




%% Forced water template for programmatic generation of MESMA setttings file
% See the example: C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\ForcedWater.txt
MesmaSettings   =   char(...
'MESMA Settings: ',...
'   ',... 
libraryFile,...% Name of input library
classNameField,...% Metadata field with end-member names
[num2str(minEmFrac) ', ' num2str(maxEmFrac) ', ' num2str(minShadeEmFrac) ', ' ...
num2str(maxShadeEmFrac) ', ' num2str(maxRMSE) ', -9999, -9999'],...% List of constraints:  min + max endmember fraction, min + max shade fraction, max rmse, residual reflectance threshold + max number of consecutive bands exceeding threshold. set value to -9999 if not used.% fullfile(dataDir,[siteDateCode 'SpecSub.tif']),...% Name of input image
fullfile(dataDir,[siteDateCode 'simSpec.tif']),...% Name of input image% fullfile(dataDir,[siteDateCode 'SpecSubMESMA']),...% Root name of MESMA output files
fullfile(dataDir,[siteDateCode 'simSpecMESMA']),...% Root name of MESMA output files
'True',...% Logical flag to open result in QGIS
'1',...% Spectral library reflectance scale factor
'0.007',...% Fusion threshold (not used)
shadeFile,...% Try name of non-photometric shade here
'1',...% Reflectance scale factor for shade end-member
'10000',...% Image reflectance scale factor
'False',...% Logical flag to create residual image
'False',...% Logical flag to use spectral weighting
'False',...% Logical flag to use band selection (Stable Zone Unmixing)
'0.99',...% Correlation threshold for Stable Zone Unmixing
'0.01',...% Correlation decrease for Stable Zone Unmixing
'x-EM Levels:',...% The rest sets up which models and end-members will be used and shouldn't change as long as the number of end-members doesn't change and we don't decide to exclude some end members on a case by case basis
'False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False',...
'classes per x-EM Level:',...
'2-EMs: True, True, True, True, True, True, True, True, True, True, True, True, True, True',...
'3-EMs: True, True, True, True, True, True, True, True, True, True, True, True, True, True',...
'4-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'5-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'6-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'7-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'8-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'9-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'10-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'11-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'12-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'13-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'14-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'15-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
'models per x-EM Level:',...
'2-EMs: True, True, True, True, True, True, True, True, True, True, True, True, True, True',...
'3-EMs: False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, False, False, False, True, False, False, False, False, False, True, False, False, False, False, True, False, False, False, True, False, False, True, False, True, True');
MesmaSettings   =   cellstr(MesmaSettings);
% writecell(MesmaSettings,'tmp.txt')
fid             =   fopen(fullfile(dataDir,[siteDateCode 'MesmaSettings.txt']),'w');
formatSpec      =   '%s\n';
for i = 1:length(MesmaSettings)
    fprintf(fid,formatSpec,MesmaSettings{i});
end
fclose(fid);
type 'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823MesmaSettings.txt'
% Looks like this worked and can be used to do a run in QGIS


%% **** 6/7/2021: PROGRAMMATIC CALL TO PYTHON COMMAND LINE IMPLEMENTATION OF MESMA
% See guidance from Stephen Hundt and his mesma setup code in:
%   C:\Users\cjl\OneDrive - DOI\HABs\SMASH\MESMA_setup
% Specifically the document Instructions for MESMA.docx
% Made a copy of his original magic.bat file and then went through the initial
% setup outlined in the instructions, which ran OK.

% But note that you have to execute magic.bat each time you want to use the
% Python command line version of mesma and then use mesma from a DOS command
% prompt initiated by magic.bat.

% So, let's try to programmatically call magic.bat and maybe add in the ability
% to issue a command to actually run the mesma

% The simplest approach is to just issue a command from MATLAB to the operating
% system using an exclamation point, like this:
    % ! C:\Users\cjl\OneDrive - DOI\HABs\SMASH\MESMA_setup\magic.bat
    
% Or we can specify a command as a string and then call it using the system
% function in MATLAB
command =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\MESMA_setup\magic.bat &';
status  =   system(command);

% Create a string to call the mesma function from the CMD window that opens up
mesmaCmd=   'mesma Library\Algae+Water.sli Class Owasco\OW20200823SpecSub.tif';
disp(mesmaCmd)
% So now we just paste this to the CMD window
% Looks like this ran OK, so now we can work on programmatically testing
% different parameter settings etc. For documentation on the MESMA python
% command line interface, see:
%   https://mesma.readthedocs.io/en/latest/userguide/mesma_cli.html

% Could even try programmatically generating a new magic.bat file that tacks on
% the command to mesma at the end to take the DOS prompt out of the loop and
% keep everything within MATLAB


%% Python mesma documentation copied from command line after running mesma -h
% usage: mesma [-h] [-l  ...]] [-r] [-f] [-u] [--min-fraction]
%              [--max-fraction] [--min-shade-fraction]
%              [--max-shade-fraction] [--max-rmse] [--residual-constraint]
%              [--residual-constraint-values  ...]] [-s] [-a] [-t]
%              [-o] [-d] [--spectral-weighing] [--band-selection]
%              [--band-selection-values  ...]] [-c]
%              library class image
% 
% Multiple Endmember Signal Mixture Analysis: calculate SMA for multiple
% endmember combinations and select the best fit based on the lowest RMSE.
% Citations: MESMA: Roberts, D.A., Gardner, M., Church, R., Ustin, S., Scheer,
% G., Green, R.O., 1998, Mapping Chaparral in the Santa Monica Mountains using
% Multiple Endmember Spectral Mixture Models, Remote Sensing of Environment, 65,
% p. 267-279. Multilevel fusion: Roberts, D.A., Dennison, P.E., Gardner, M.,
% Hetzel, Y., Ustin, S.L., Lee, C., 2003, Evaluation of the Potential of
% Hyperion for Fire Danger Assessment by Comparison to the Airborne
% Visible/Infrared Imaging Spectrometer, IEEE Transactions on Geoscience and
% Remote Sensing, 41, p. 1297-1310. Spectral band weighing: Somers, B.,
% Delalieux, S, Stuckens, J , Verstraeten, W.W, Coppin, P., 2009, A weighted
% linear spectral mixture analysis approach to address endmember variability in
% agricultural production systems, International Journal of Remote Sensing, 30,
% p. 139-147. Spectral band selection: Somers, B., Delalieux, S., Verstraeten,
% W.W., van Aardt, J.A.N., Albrigo, G., Coppin, P., 2010, An automated waveband
% selection technique for optimized hyperspectral mixture analysis.
% International Journal of Remote Sensing, 31, p. 5549-5568.
% 
% positional arguments:
%   library               spectral library file
%   class                 metadata header with the spectrum classes
%   image                 image file (singular)
% 
% optional arguments:
%   -h, --help            show this help message and exit
%   -l  ...], --complexity-level  ...]
%                         the complexity levels for unmixing. e.g. 2 3 4 for 2-,
%                         3- and 4-EM models (default: 2 3)
%   -r, --reflectance-scale-library
%                         library reflectance scale factor (default: derived
%                         from data as 1, 1 000 or 10 000)
%   -f, --fusion-threshold
%                         Models with a higher number of classes (e.g. 4-EM vs.
%                         3-EM models) are only chosen when the decrease in RMSE
%                         is larger than this threshold (default: 0.007)
%   -u, --unconstrained   run mesma without constraints (default off)
%   --min-fraction      minimum allowable endmember fraction (default -0.05),
%                         use -9999 to set no constraint
%   --max-fraction      maximum allowable endmember fraction (default 1.05),
%                         use -9999 to set no constraint
%   --min-shade-fraction
%                         minimum allowable shade fraction (default 0.00), use
%                         -9999 to set no constraint
%   --max-shade-fraction
%                         maximum allowable shade fraction (default 0.80), use
%                         -9999 to set no constraint
%   --max-rmse          maximum allowable RMSE (default 0.025), use -9999 to
%                         set no constraint
%   --residual-constraint
%                         use a residual constraint (default off)
%   --residual-constraint-value  ...]
%                         two values (residual threshold, number of bands): the
%                         number of consecutive bands that the residual values
%                         are allowed to exceed the given threshold (default:
%                         0.025 7)
%   -s, --reflectance-scale-image
%                         image reflectance scale factor (default: derived from
%                         data as 1, 1 000 or 10 000)
%   -a, --shade       non-photometric shade, saved as library file with one
%                         spectrum
%   -t, --reflectance-scale-shade
%                         shade reflectance scale factor (default: derived from
%                         data as 1, 1 000 or 10 000)
%   -o, --output      output image file (default: in same folder as image
%                         with extension '_mesma_datetime'
%   -d, --residuals-image
%                         create a residual image as output (default off)
%   --spectral-weighing   use spectral weighing algorithm (default off)
%   --band-selection      use band selection algorithm (default off)
%   --band-selection-values  ...]
%                         two values: correlation threshold and correlation
%                         decrease (default: 0.99 0.01)
%   -c, --cores       the number of logical cpu cores available for this
%                         process (default: 1)


%% Here's the original magic.bat file
% ::QGIS installation folder
% set OSGEO4W_ROOT="C:\OSGeo4W64"
% 
% ::set defaults, clean path, load OSGeo4W modules (incrementally)
% call %OSGEO4W_ROOT%\bin\o4w_env.bat
% call qt5_env.bat
% call py3_env.bat
% 
% ::lines taken from python-qgis.bat
% set QGIS_PREFIX_PATH=%OSGEO4W_ROOT%\apps\qgis
% set PATH=%QGIS_PREFIX_PATH%\bin;%PATH%
% 
% ::make PyQGIS packages available to Python
% set PYTHONPATH=%QGIS_PREFIX_PATH%\python;%PYTHONPATH%
% 
% :: GDAL Configuration (https://trac.osgeo.org/gdal/wiki/ConfigOptions)
% :: Set VSI cache to be used as buffer, see #6448 and
% set GDAL_FILENAME_IS_UTF8=YES
% set VSI_CACHE=TRUE
% set VSI_CACHE_SIZE=1000000
% set QT_PLUGIN_PATH=%QGIS_PREFIX_PATH%\qtplugins;%OSGEO4W_ROOT%\apps\qt5\plugins
% 
% ::enable/disable QGIS debug messages
% set QGIS_DEBUG=1
% 
% ::open the OSGeo4W Shell
% ::@echo on
% ::@if [%1]==[] (echo run o-help for a list of available commands & cmd.exe /k) else (cmd /c "%*")
% mesma -h
% ::Should be able to add our customized command from MATLAB here to actually run the mesma, then end with exit
% mesma Library\Algae+Water.sli Class Owasco\OW20200823SpecSub.tif
% exit


%% Set parameters we will pass to the MESMA core algorithm
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\Algae+Water.sli';
% Metadata field with end-member names
classNameField  =   'Class';
% Name of image file
imgSubFile      =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823specSub.tif';
% Scale factors for library and image
libScaleFactor  =   1;
imgScaleFactor  =   10000;
% Constraints
%   Minimum end member fraction
minEmFrac       =   -0.05;
%   Maximum end member fraction
maxEmFrac       =   1.05;
%   Minimum shade end member fraction
minShadeEmFrac  =   -0.01;
%   Maximum shade end member fraction
maxShadeEmFrac  =   0.01;
%   Maximum RMSE
maxRMSE         =   0.025;
% Name of shade end member file
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\DetroitLakeWaterOnly.sli';


%% TO DO: Work on lookup table and class per em to force water
% Set up look_up_table and em_per_class inputs to specify which models and EM's to use
% From the mesma.py documentation:
    % look_up_table – all endmember combinations (=models) for MESMA; ordered
    %                 per complexity level and per class-model; n_models x n_endmembers
    % em_per_class – a list of all library indices per endmember class

% Example from the QGIS input file we used previously to force water models    
% 'x-EM Levels:',...% The rest sets up which models and end-members will be used and shouldn't change as long as the number of end-members doesn't change and we don't decide to exclude some end members on a case by case basis
% 'False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False',...
% 'classes per x-EM Level:',...
% '2-EMs: True, True, True, True, True, True, True, True, True, True, True, True, True, True',...
% '3-EMs: True, True, True, True, True, True, True, True, True, True, True, True, True, True',...
% '4-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% '5-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% '6-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% '7-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% '8-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% '9-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% '10-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% '11-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% '12-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% '13-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% '14-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% '15-EMs: False, False, False, False, False, False, False, False, False, False, False, False, False, False',...
% 'models per x-EM Level:',...
% '2-EMs: True, True, True, True, True, True, True, True, True, True, True, True, True, True',...
% '3-EMs: False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, False, False, False, True, False, False, False, False, False, True, False, False, False, False, True, False, False, False, True, False, False, True, False, True, True');

% For the library with 13 algal taxa and water (14 EM), the table will be 14 EM
% X 15 models, so for 2 and 3-EM models we have
lut     =   char(...
    'True, True, True, True, True, True, True, True, True, True, True, True, True, True',...%2-EM models
    'True, True, True, True, True, True, True, True, True, True, True, True, True, True',...%3-EM models
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False',...%4-EMs
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False',...%5-EMs 
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False',...%6-EMs
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False',...%7-EMs
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False',...%8-EMs: 
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False',...%9-EMs
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False',...%10-EMs
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False',...%11-EMs
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False',...%12-EMs
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False',...%13-EMs
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False',...%14-EMs
    'False, False, False, False, False, False, False, False, False, False, False, False, False, False');%15-EMs

% em_per_class – a list of all library indices per endmember class
% For the three-end member models forced to include water
% '3-EMs: False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, True, False, False, False, False, False, False, True, False, False, False, False, False, True, False, False, False, False, True, False, False, False, True, False, False, True, False, True, True');


%% Programmatically generate new magic.bat file with actual call to mesma embedded
% Set up the actual call to the mesma routine itself first so we can embed it in
% the .bat file we're writing here
% Create a unique file name for the .bat file to be generated
timeStamp   =   datestr(now,'yyyymmdd_HHMMSS');
batFileName =   fullfile(dataDir,[siteDateCode '_CallMesma_' timeStamp '.bat']);
outputRoot  =   fullfile(dataDir,[siteDateCode '_MesmaOutput_' timeStamp]);
mesmaCmd    =   ['mesma "' libraryFile '" ' classNameField ' "' imgSubFile '"' ...
                 ' -l 2' ... Set complexity level to run all 2-EM models
                 ' -f 0.007' ... The default fusion threshold is 0.007. Models with a higher number of classes (e.g. 4-EM vs. 3-EM models) are only chosen when the decrease in RMSE is larger than this threshold
                 ' --min-fraction ' num2str(minEmFrac) .... By default, -0.05 (minimum endmember fraction)
                 ' --max-fraction ' num2str(maxEmFrac) .... By default, 1.05 (maximum endmember fraction)
                 ' --min-shade-fraction ' num2str(minShadeEmFrac) .... By default, 0.00 (minimum shade fraction)
                 ' --max-shade-fraction ' num2str(maxShadeEmFrac) .... By default, 0.80 (maximum shade fraction)
                 ' --max-rmse ' num2str(maxRMSE) .... By default, 0.025 (maximum RMSE)
                 ' -r ' num2str(libScaleFactor) ... reflectance scale factor of the spectral library
                 ' -s ' num2str(imgScaleFactor) ... reflectance scale factor of the image
                 ' -o "' outputRoot '"']; % root file name for all the output
                 
% Now for the preamble to setup mesma and call it from the DOS command line
CallMesma       =   char(...
    '::QGIS installation folder',...
    'set OSGEO4W_ROOT="C:\OSGeo4W64"',...
    '   ',...
    '::set defaults, clean path, load OSGeo4W modules (incrementally)',...
    'call %OSGEO4W_ROOT%\bin\o4w_env.bat',...
    'call qt5_env.bat',...
    'call py3_env.bat',...
    '   ',...
    '::lines taken from python-qgis.bat',...
    'set QGIS_PREFIX_PATH=%OSGEO4W_ROOT%\apps\qgis',...
    'set PATH=%QGIS_PREFIX_PATH%\bin;%PATH%',...
    '   ',...
    '::make PyQGIS packages available to Python',...
    'set PYTHONPATH=%QGIS_PREFIX_PATH%\python;%PYTHONPATH%',...
    '   ',...
    ':: GDAL Configuration (https://trac.osgeo.org/gdal/wiki/ConfigOptions)',...
    ':: Set VSI cache to be used as buffer, see #6448 and',...
    'set GDAL_FILENAME_IS_UTF8=YES',...
    'set VSI_CACHE=TRUE',...
    'set VSI_CACHE_SIZE=1000000',...
    'set QT_PLUGIN_PATH=%QGIS_PREFIX_PATH%\qtplugins;%OSGEO4W_ROOT%\apps\qt5\plugins',...
    '   ',...
    '::enable/disable QGIS debug messages',...
    'set QGIS_DEBUG=1',...
    '   ',...
    '::open the OSGeo4W Shell',...
    '   ',...
    mesmaCmd,...
    'exit');
CallMesma       =   cellstr(CallMesma);
fid             =   fopen(batFileName,'w');
formatSpec      =   '%s\n';
for i = 1:length(CallMesma)
    fprintf(fid,formatSpec,CallMesma{i});
end
fclose(fid);
% type 'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823CallMesma.bat'


%% Now call the new .bat file to run the mesma
% Specify the command as a string and then call it using the system function
command =   [batFileName ' &'];
status  =   system(command);
if status == 0
    disp('MESMA executed successfully, now importing results ...')
else
    disp('WARNING: MESMA did not execute as intended ...')
end


%% Next move on to programmatically importing the MESMA results


%% Use Hyperspectral Image Processing Toolbox to import hypercube and display RGB
Hcube       =   hypercube(fullfile(dataDir,[siteDateCode 'SpecSub.tif']),wvlSub);
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


%% Get list of taxa from library and display numeric codes
Taxa    =   readtable([libraryFile(1:end-3) 'csv']);
Taxa    =   Taxa.Class;
TaxaNums=   cell(size(Taxa));
for i = 1:length(Taxa)
    TaxaNums{i} =   [num2str(i) ': ' Taxa{i}];
end
disp('Algal taxa numeric codes and taxa names:')
disp(TaxaNums)
hFig    =   uifigure('Position',[50 600 250 300],'Name','Algal Taxa Listing');
hLabel  =   uilistbox(hFig,'Items',TaxaNums,'Position',[5 5 240 290]);


%% *MESMA file contains the best model for each pixel as library spectra number per class
% This file has the best model [nb of bands = nb of classes] - each band
% contains the library spectra number per class
% Value of unmodeled pixels in output: -1
% Value of pixels with no data in output: -2
% outputRoot  =   fullfile(dataDir,[siteDateCode '_MesmaOutput_' timeStamp]);
% modelFile       =   fullfile(dataDir,[siteDateCode 'SpecSubMESMA']);
modelFile       =   fullfile(dataDir,[siteDateCode '_MesmaOutput_' timeStamp]);
model           =   envi2matlab(modelFile);
% Flip to account for geo-referencing
model           =   flipud(model);
% Reset nodata pixels to NaN
model(model==-2)=   NaN;
% % Reset unmodeled pixels to NaN;
% model(model==-1)=   NaN;
% Apply mask
model           =   applyMask(model,waterMask);


%% Display classified map with a handmade legend - might not be right for 3-EM models
modelSum    =   zeros(size(model,[1 2]));
for i = 1:size(model,3)
    iTaxa       =   find(model(:,:,i)==i-1);
    tmp         =   zeros(size(model,[1 2]));
    tmp(iTaxa)  =   i; %#ok<FNDSB>
    modelSum    =   modelSum + tmp;
end
modelSum(~waterMask)    =   NaN;
figure
imagesc(modelSum,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask | modelSum==0);
map             =   colormap(lines(size(model,3)+1));
hCbar           =   colorbar;
hCbar.Ticks     =   [];
% A couple of hard-wired parameters for the legend layout
legLabelOffset  =   -0.1;
legBoxWidth     =   0.05;
legBoxEdge      =   hCbar.Position(1)+legLabelOffset;
legBoxHeight    =   hCbar.Position(4)/(length(Taxa)+1);
legBoxBot       =   hCbar.Position(2);
Legend          =   [{'Unclassified'}; Taxa];
for i = 1:length(Legend)
    textBoxDim  =   [legBoxEdge legBoxBot+legBoxHeight*(i-1) ...
                     legBoxWidth legBoxHeight];
    annotation('textbox',textBoxDim,'String',Legend(i),'FitBoxToText','on',...
               'fontname','arial','fontsize',12,'linestyle','none',...
               'verticalalignment','bottom','horizontalalignment','left');
end
axis equal; axis xy
xlabel('Easting (m)')
ylabel('Northing (m)')
title([siteDateCode ': MESMA-based classification']);


%% Specify taxa number to display area mapped as a single class
taxaId          =   input('Specify taxa ID number to display locations classified as that taxa: ');
iTaxa           =   find(model(:,:,taxaId)==taxaId-1);
classMask       =   false(size(model(:,:,taxaId)));
classMask(iTaxa)=   true;
figure
imagesc(model(:,:,taxaId),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & classMask);
crameri('-bamako')    
axis image; axis xy; axis on;
xlabel('Easting (m)')
ylabel('Northing (m)')
title(['MESMA pixels with best model for taxa: ' Taxa(taxaId)])


%% *MESMA_fractions file contains the model's fractions, including shade
% The model’s fractions [nb of bands = nb of classes + 1], including a shade fraction
% Value of unmodeled pixels in output: 0
% Value of pixels with no data in output: 0
% fractionsFile   =   fullfile(dataDir,[siteDateCode 'SpecSubMESMA_fractions']);
fractionsFile   =   fullfile(dataDir,[siteDateCode '_MesmaOutput_' timeStamp '_fractions']);
fractions       =   envi2matlab(fractionsFile);
% Flip to account for geo-referencing
fractions       =   flipud(fractions);
% Reset nodata pixels to NaN
fractions(fractions==0)=   NaN;
% Reset unmodeled pixels to NaN;
fractions(fractions==0)=   NaN;
% Apply mask
fractions       =   applyMask(fractions,waterMask);


%% Specify taxa number to display and display fraction image for that taxa
taxaId          =   input('Specify taxa ID number to display fractions for that taxa: ');
figure
imagesc(fractions(:,:,taxaId),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & ~isnan(fractions(:,:,taxaId)));
axis image; axis xy; axis on; 
crameri('-bamako');
colorbar
xlabel('Easting (m)')
ylabel('Northing (m)')
title(['MESMA fractions for taxa: ' Taxa(taxaId)])


%% *MESMA_rmse file contains the model's RMSE
% The model’s RMSE [nb of bands = 1]
% Value of unmodeled pixels in output: 9999
% Value of pixels with no data in output: 9998
% rmseFile   =   fullfile(dataDir,[siteDateCode 'SpecSubMESMA_rmse']);
rmseFile    =   fullfile(dataDir,[siteDateCode '_MesmaOutput_' timeStamp '_rmse']);
rmse       =   envi2matlab(rmseFile);
% Flip to account for geo-referencing
rmse       =   flipud(rmse);
% Reset nodata pixels to NaN
rmse(rmse==9998)=   NaN;
% Reset unmodeled pixels to NaN;
rmse(rmse==9999)=   NaN;
% Apply mask
rmse       =   applyMask(rmse,waterMask);


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







%% **** 6/18/2021: REVISIT OWASCO WITH WATER END MEMBER AS SHADE
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
load SMASHowasco.mat


%% Set parameters we will pass to the MESMA core algorithm
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\AlgaeTaxaOnly.sli';
% Metadata field with end-member names
classNameField  =   'Class';
% Name of image file
imgSubFile      =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823specSub.tif';
% Scale factors for library and image
libScaleFactor  =   1;
imgScaleFactor  =   10000;
% Constraints
%   Minimum end member fraction
minEmFrac       =   -0.05;
%   Maximum end member fraction
maxEmFrac       =   1.05;
%   Minimum shade end member fraction
minShadeEmFrac  =   -0.05;
%   Maximum shade end member fraction
maxShadeEmFrac  =   1.05;
%   Maximum RMSE
maxRMSE         =   0.025;
% Name of shade end member file
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\DetroitLakeWaterOnly.sli';


%% Programmatically generate new magic.bat file with actual call to mesma embedded
% Set up the actual call to the mesma routine itself first so we can embed it in
% the .bat file we're writing here
% Create a unique file name for the .bat file to be generated
timeStamp   =   datestr(now,'yyyymmdd_HHMMSS');
batFileName =   fullfile(dataDir,[siteDateCode '_CallMesma_' timeStamp '.bat']);
outputRoot  =   fullfile(dataDir,[siteDateCode '_MesmaOutput_' timeStamp]);
mesmaCmd    =   ['mesma "' libraryFile '" ' classNameField ' "' imgSubFile '"' ...
                 ' -l 2 3' ... Set complexity level to run all 2 and 3-EM models
                 ' -f 0.007' ... The default fusion threshold is 0.007. Models with a higher number of classes (e.g. 4-EM vs. 3-EM models) are only chosen when the decrease in RMSE is larger than this threshold
                 ' --min-fraction ' num2str(minEmFrac) .... By default, -0.05 (minimum endmember fraction)
                 ' --max-fraction ' num2str(maxEmFrac) .... By default, 1.05 (maximum endmember fraction)
                 ' --min-shade-fraction ' num2str(minShadeEmFrac) .... By default, 0.00 (minimum shade fraction)
                 ' --max-shade-fraction ' num2str(maxShadeEmFrac) .... By default, 0.80 (maximum shade fraction)
                 ' --max-rmse ' num2str(maxRMSE) .... By default, 0.025 (maximum RMSE)
                 ' -a "' shadeFile '"' ... shade file
                 ' -r ' num2str(libScaleFactor) ... reflectance scale factor of the spectral library
                 ' -s ' num2str(imgScaleFactor) ... reflectance scale factor of the image
                 ' -o "' outputRoot '"']; % root file name for all the output
                 
% Now for the preamble to setup mesma and call it from the DOS command line
CallMesma       =   char(...
    '::QGIS installation folder',...
    'set OSGEO4W_ROOT="C:\OSGeo4W64"',...
    '   ',...
    '::set defaults, clean path, load OSGeo4W modules (incrementally)',...
    'call %OSGEO4W_ROOT%\bin\o4w_env.bat',...
    'call qt5_env.bat',...
    'call py3_env.bat',...
    '   ',...
    '::lines taken from python-qgis.bat',...
    'set QGIS_PREFIX_PATH=%OSGEO4W_ROOT%\apps\qgis',...
    'set PATH=%QGIS_PREFIX_PATH%\bin;%PATH%',...
    '   ',...
    '::make PyQGIS packages available to Python',...
    'set PYTHONPATH=%QGIS_PREFIX_PATH%\python;%PYTHONPATH%',...
    '   ',...
    ':: GDAL Configuration (https://trac.osgeo.org/gdal/wiki/ConfigOptions)',...
    ':: Set VSI cache to be used as buffer, see #6448 and',...
    'set GDAL_FILENAME_IS_UTF8=YES',...
    'set VSI_CACHE=TRUE',...
    'set VSI_CACHE_SIZE=1000000',...
    'set QT_PLUGIN_PATH=%QGIS_PREFIX_PATH%\qtplugins;%OSGEO4W_ROOT%\apps\qt5\plugins',...
    '   ',...
    '::enable/disable QGIS debug messages',...
    'set QGIS_DEBUG=1',...
    '   ',...
    '::open the OSGeo4W Shell',...
    '   ',...
    mesmaCmd,...
    'exit');
CallMesma       =   cellstr(CallMesma);
fid             =   fopen(batFileName,'w');
formatSpec      =   '%s\n';
for i = 1:length(CallMesma)
    fprintf(fid,formatSpec,CallMesma{i});
end
fclose(fid);
% type 'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823CallMesma.bat'


%% Now call the new .bat file to run the mesma
% Specify the command as a string and then call it using the system function
command =   [batFileName ' &'];
status  =   system(command);
if status == 0
    disp('MESMA executed successfully, now importing results ...')
else
    disp('WARNING: MESMA did not execute as intended ...')
end


%% Next move on to programmatically importing the MESMA results


%% Use Hyperspectral Image Processing Toolbox to import hypercube and display RGB
Hcube       =   hypercube(fullfile(dataDir,[siteDateCode 'SpecSub.tif']),wvlSub);
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


%% Get list of taxa from library and display numeric codes
Taxa    =   readtable([libraryFile(1:end-3) 'csv']);
Taxa    =   Taxa.Class;
TaxaNums=   cell(size(Taxa));
for i = 1:length(Taxa)
    TaxaNums{i} =   [num2str(i) ': ' Taxa{i}];
end
disp('Algal taxa numeric codes and taxa names:')
disp(TaxaNums)
hFig    =   uifigure('Position',[50 600 250 300],'Name','Algal Taxa Listing');
hLabel  =   uilistbox(hFig,'Items',TaxaNums,'Position',[5 5 240 290]);


%% *MESMA file contains the best model for each pixel as library spectra number per class
% This file has the best model [nb of bands = nb of classes] - each band
% contains the library spectra number per class
% Value of unmodeled pixels in output: -1
% Value of pixels with no data in output: -2
modelFile       =   fullfile(dataDir,[siteDateCode '_MesmaOutput_' timeStamp]);
model           =   envi2matlab(modelFile);
% Flip to account for geo-referencing
model           =   flipud(model);
% Reset nodata pixels to NaN
model(model==-2)=   NaN;
% Reset unmodeled pixels to NaN;
model(model==-1)=   NaN;
% Apply mask
model           =   applyMask(model,waterMask);


%% Display classified map with a handmade legend - might not be right for 3-EM models
modelSum    =   zeros(size(model,[1 2]));
for i = 1:size(model,3)
    iTaxa       =   find(model(:,:,i)==i-1);
    tmp         =   zeros(size(model,[1 2]));
    tmp(iTaxa)  =   i; %#ok<FNDSB>
    modelSum    =   modelSum + tmp;
end
modelSum(~waterMask)    =   NaN;
figure
imagesc(modelSum,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask | modelSum==0);
map             =   colormap(lines(size(model,3)+1));
hCbar           =   colorbar;
hCbar.Ticks     =   [];
% A couple of hard-wired parameters for the legend layout
legLabelOffset  =   -0.1;
legBoxWidth     =   0.05;
legBoxEdge      =   hCbar.Position(1)+legLabelOffset;
legBoxHeight    =   hCbar.Position(4)/(length(Taxa)+1);
legBoxBot       =   hCbar.Position(2);
Legend          =   [{'Unclassified'}; Taxa];
for i = 1:length(Legend)
    textBoxDim  =   [legBoxEdge legBoxBot+legBoxHeight*(i-1) ...
                     legBoxWidth legBoxHeight];
    annotation('textbox',textBoxDim,'String',Legend(i),'FitBoxToText','on',...
               'fontname','arial','fontsize',12,'linestyle','none',...
               'verticalalignment','bottom','horizontalalignment','left');
end
axis equal; axis xy
xlabel('Easting (m)')
ylabel('Northing (m)')
title([siteDateCode ': MESMA-based classification']);


%% Specify taxa number to display area mapped as a single class
taxaId          =   input('Specify taxa ID number to display locations classified as that taxa: ');
iTaxa           =   find(model(:,:,taxaId)==taxaId-1);
classMask       =   false(size(model(:,:,taxaId)));
classMask(iTaxa)=   true;
figure
imagesc(model(:,:,taxaId),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & classMask);
crameri('-bamako');    
axis image; axis xy; axis on;
xlabel('Easting (m)')
ylabel('Northing (m)')
title(['MESMA pixels with best model for taxa: ' Taxa(taxaId)])


%% *MESMA_fractions file contains the model's fractions, including shade
% The model’s fractions [nb of bands = nb of classes + 1], including a shade fraction
% Value of unmodeled pixels in output: 0
% Value of pixels with no data in output: 0
% fractionsFile   =   fullfile(dataDir,[siteDateCode 'SpecSubMESMA_fractions']);
fractionsFile   =   fullfile(dataDir,[siteDateCode '_MesmaOutput_' timeStamp '_fractions']);
fractions       =   envi2matlab(fractionsFile);
% Flip to account for geo-referencing
fractions       =   flipud(fractions);
% Reset nodata pixels to NaN
fractions(fractions==0)=   NaN;
% Reset unmodeled pixels to NaN;
fractions(fractions==0)=   NaN;
% Apply mask
fractions       =   applyMask(fractions,waterMask);


%% Specify taxa number to display and display fraction image for that taxa
taxaId          =   input('Specify taxa ID number to display fractions for that taxa: ');
figure
imagesc(fractions(:,:,taxaId),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & ~isnan(fractions(:,:,taxaId)));
axis image; axis xy; axis on; 
crameri('-bamako');
colorbar
xlabel('Easting (m)')
ylabel('Northing (m)')
title(['MESMA fractions for taxa: ' Taxa(taxaId)])


%% Looking at all 14 end members (including water as shade) shows mostly micro, some merismo, so make RGB
taxaId1         =   input('Specify taxa ID number for fraction to display as red: ');
taxaId2         =   input('Specify taxa ID number for fraction to display as green: ');
% Assume water (shade) is the last end member in the fractions stack
rgbFrac         =   fractions(:,:,[taxaId1 taxaId2 size(fractions,3)]);
figure
% imshow(rgbFrac)
imagesc(rgbFrac,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',waterMask & isnan);
axis image; axis xy; axis on; 
xlabel('Easting (m)')
ylabel('Northing (m)')
title({'Fraction composite:',['R:' Taxa{taxaId1} ' G:' Taxa{taxaId2} ...
       ' B:Water']})


%% *MESMA_rmse file contains the model's RMSE
% The model’s RMSE [nb of bands = 1]
% Value of unmodeled pixels in output: 9999
% Value of pixels with no data in output: 9998
% rmseFile   =   fullfile(dataDir,[siteDateCode 'SpecSubMESMA_rmse']);
rmseFile    =   fullfile(dataDir,[siteDateCode '_MesmaOutput_' timeStamp '_rmse']);
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


%% Save what we have so far
clear ans fid h* i iTaxa MesmaSettings status tmp
save SMASHowasco.mat