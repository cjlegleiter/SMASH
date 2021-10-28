function [R,ndwi,waterMask,ndci,ci,wvlSub,imgSpecFiltSmooth,imgSubFile] = ...
    preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput)
% Prepare inputs for MESMA of algal taxa as part of HAB remote sensing workflow
%
%% preSMASH.m:
%   Prepare inputs for MESMA of algal taxa as the first steps in a workflow for
%   remote sensing of HABs.  This function imports the original hyperspectral
%   image, prompts the user to select a spatial subset (region of interest, or
%   ROI), calculates the Normalized Difference Water Index (NDWI), prompts the
%   user to interactively create a water mask (or an existing mask as a
%   shapefile can be used instead), calculates the Normalized Difference
%   Chlorophyll Index (NDCI), subsets the image spectrally to the desired
%   wavelength range, applies a Savitzky-Golay filter to the spectrum for each
%   pixel and a Weiner spatial smoothing filter to the image, and exports then
%   exports the final, pre-processed image for further processing (i.e., MESMA).
%   The outputs from this function are used as inputs to runSMASH.m.
%
%% SYNTAX:
% [R,ndwi,waterMask,ndci,ci,wvlSub,imgSpecFiltSmooth,imgSubFile] = ...
%     preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput);
%
%% INPUTS:
%   dataDir:    String specifying directory where outputs from this function
%               will be created
%   siteCode:   String with two-letter site code
%   dateCode:   String with 8-digit date code (YYYYMMDD)
%   epsgCode:   Scalar with EPSG projection code for writing GEO-TIFF output;
%               see the EPSG website at https://epsg.org/home.html
%   wvl2keep:   2-element vector with wavelength range to keep in subset: 
%               [minWvl maxWvl] (units assumed to be nm)
%   imgFile:    String with full path to original image file
%   maskInput:  Use an empty array (i.e., []) to interactively define a mask
%               based on the NDWI image, or provide a string with the full path
%               to an existing water mask represented as a shapefile
% 
%% OUTPUTS:
%   R:          Geo-referencing object
%   ndwi:       Normalized Difference Water Index (NDWI) output image
%   waterMask:  Binary water mask image
%   ndci:       Normalized Difference Chlorophyll Index (NDCI) output image
%   ci:         Cyanobatrerial index (CI) output image
%   wvlSub:     Vector of wavelengths retained (units of nm)
%   imgSpecFiltSmooth: Spatially and spectrally subset image with Savitzky-Golay
%               spectral filter and spatial smoothing filter applied
%   imgSubFile: String with file name of output image to use in MESMA
%
%% NOTES:
% > Requires the MATLAB Image Processing and Mapping Toolboxes
% > The Hyperspectral Imaging Toolbox add-in must also be installed
% > Code contains comments with inputs used for Lake Owasco prototype
% > The Savitzky-Golay spectral filter parameters are hardwired in the code to
%   a 3rd-order polynomial with a window size of 7, applied twice
% > The window size of the Wiener spatial smoothing filter is hardwired as 3 X 3
%
%% FUNCTION SUMMARY:
% [R,ndwi,waterMask,ndci,wvlSub,imgSpecFiltSmooth,imgSubFile] = ...
%     preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput);

%% CREDITS:
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 06/18/2021
% 08/12/2021 - Updated to include spectral filtering and spatial smoothing
% 08/18/2021 - Updated to calculate cyanobacterial index (CI) image
% Also see prototyping code in: C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SMASHowasco.m
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\preSMASH.m


%% Set data directory, two-letter site code, date, EPSG projection, and bands to keep
% cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
% dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco';
% siteCode    =   'OW';
% dateCode    =   '20200823';
% % WGS 84 UTM Zone 18N is EPSG 32618;
% epsgCode    =   32618;
siteDateCode=   [siteCode dateCode];
% % Specify bands to keep to match convolved library
% iKeepWvl    =   1:158;


%% Import raw image file, get geo-referencing info, and display a NIR band
% imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\rawImage\DESIS-HSI-L2A-DT0489826624_003-20200823T165204-V0210-SPECTRAL_IMAGE.tif';
hdrFile =   [imgFile(1:end-3) 'hdr'];
filepath=   fileparts(imgFile);
wvlFile =   fullfile(filepath,'wvlList.txt');
if exist(hdrFile,'file')
    % Read the header, mainly to get the wavelengths
    [~,~,~,~,~,~,~,~,wvl]   =   readEnviHeader(imgFile);
elseif exist(wvlFile,'file')
    % Or if there's no header, use the wavelength list text file we made for DESIS
    wvl                     =   load(wvlFile);
end
% Use a function from the hyperspectral imaging toolbox to read the full data cube            
hcube                   =   hypercube(imgFile,wvl);
% Flip to account for geo-referencing
img                     =   flipud(hcube.DataCube);
clear hcube

% Read a single band with readgeoraster to get spatial referencing information
[~,iNirBand]            =   min(abs(wvl-750));
[nirBand,R]             =   readgeoraster(imgFile,'Bands',iNirBand);
% Flip to account for geo-referencing
nirBand                 =   flipud(nirBand);
% Display with a 98% contrast stretch
figure
imshow(nirBand,prctile(nirBand(nirBand>0),[1 99]),...
    'XData',R.XWorldLimits,'YData',R.YWorldLimits);
axis xy; axis on
xlabel('Easting (m)')
ylabel('Northing (m)')
title({['Raw image for site/date: ' siteDateCode],...
       ['Displaying NIR band: ' num2str(wvl(iNirBand)) ' nm']})

   
%% Prompt user to draw rectangular region of interest (ROI)
disp('Click and drag to create a rectangular Region of Interest (ROI) ...')
roi             =   drawrectangle('Label','ROI','Color',[0 0 1]);
% Use this to get the coordinates of the ROI bounding box
xLimits         =   [roi.Position(1) roi.Position(1)+roi.Position(3)];
yLimits         =   [roi.Position(2) roi.Position(2)+roi.Position(4)];


%% Crop the image to this ROI
% Note the flipping back and forth required to make the geo-referencing work
[imgSub,Rsub]   =   mapcrop(flipud(img),R,xLimits,yLimits);
imgSub          =   flipud(imgSub);
tmp             =   imgSub(:,:,iNirBand);
% Display the cropped image
figure
imshow(tmp,prctile(tmp(tmp>0),[1 99]),'XData',Rsub.XWorldLimits,...
       'YData',Rsub.YWorldLimits);
axis xy; axis on
xlabel('Easting (m)')
ylabel('Northing (m)')
title({['ROI subset image for site/date: ' siteDateCode],...
       ['Displaying NIR band: ' num2str(wvl(iNirBand)) ' nm']})
% Overwrite original image and geo-referencing to avoid confusion and save space
img             =   imgSub;
R               =   Rsub;
clear imgSub Rsub tmp nirBand


%% Calculate NDWI and display NDWI image
% (Band 63 - band 182) / (band 63 + band 182), where 63 is green (560.6 nm) and
% 182 is NIR (865.4 nm)
[~,iGrnWvl]     =   min(abs(wvl-560.6));
[~,iNirWvl]     =   min(abs(wvl-865.4));
img             =   double(img);
ndwi            =   (img(:,:,iGrnWvl) - img(:,:,iNirWvl))./(img(:,:,iGrnWvl) + img(:,:,iNirWvl));
figure
imshow(ndwi,'XData',R.XWorldLimits,'YData',R.YWorldLimits,'DisplayRange',[]);
axis xy; axis on; colorbar
title([siteDateCode ': Normalized Difference Water Index'])
xlabel('Easting (m)')
ylabel('Northing (m)')


%% Create water mask based on NDWI (or import existing mask as shapefile) and apply to image
% maskInput   =   [];
if isempty(maskInput)
    % Note that we have to multiply by 1000 to get a reasonable data range
    waterMask   = buildWaterMask1band(ndwi*1000,R.XWorldLimits,R.YWorldLimits,1);
    % Used [-500 1000]
else
    % Allow for input of an external mask as a shapefile
    waterMask   =   vector2mask(ndwi,R.XWorldLimits,R.YWorldLimits);
end
close all

% Exclude pixels with zeros in NDWI image from the mask as well
zeroMask    =   ndwi == 0;
waterMask   =   waterMask & ~zeroMask;

% Apply mask to image
imgMask =   applyMask(img,waterMask);


%% Display the masked image for a green band
figure
imagesc(imgMask(:,:,iGrnWvl),'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
    'AlphaData',waterMask);
axis image; axis xy; axis on; colormap gray
colorbar
title([siteDateCode ': Masked image for ' num2str(wvl(iGrnWvl)) ' nm band'])
xlabel('Easting (m)')
ylabel('Northing (m)')


%% Calculate NDCI and display result
% (Band 120 - Band 104) / (Band 120 + Band 104)
% Where 120 is NIR (706.7 nm) and 104 is Red (665.3 nm)
[~,iRedWvl] =   min(abs(wvl-665.3));
[~,iNirWvl2]=   min(abs(wvl-706.7));
ndci    =   (imgMask(:,:,iNirWvl2) - imgMask(:,:,iRedWvl))./...
            (imgMask(:,:,iNirWvl2) + imgMask(:,:,iRedWvl));
ndci    =   applyMask(ndci,waterMask);
figure
imagesc(ndci,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
    'AlphaData',waterMask);
axis image; axis xy; axis on;
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


%% Spectral subset to match convolved library
iKeepWvl    =   wvl>wvl2keep(1) & wvl<wvl2keep(2);
imgSpecSub  =   img(:,:,iKeepWvl);
imgSpecSub  =   applyMask(imgSpecSub,waterMask);
wvlSub      =   wvl(iKeepWvl);


%% Apply Savitzky-Golay filter to spectrum for each pixel
% Apply smoothing filter to spectra
% Parameters for the Savitzky-Golay smoothing filter to be applied to raw data:
nFilt   =   2;
order   =   3;
window  =   7;
% Now apply this approach to the image in a double for loop over the pixels
disp('Applying a Savtzky-Golay filter to the spectrum for each pixel ...')
imgSpecFilt =   zeros(size(imgSpecSub,1),size(imgSpecSub,2),size(imgSpecSub,3));
tmp         =   imgSpecSub;
for i = 1:size(tmp,1)
    for j = 1:size(tmp,2)
        tmpPix  =   squeeze(tmp(i,j,:));
        for iFilt=1:nFilt
            tmpPix  =   sgolayfilt(tmpPix,order,window);
        end    
        imgSpecFilt(i,j,:)   =   tmpPix;
    end
end


%% Now apply a spatial smoothing filter as well
% doc wiener2
filtSize            =   [3 3];
imgSpecFiltSmooth   =   zeros(size(imgSpecFilt));
for i = 1:size(imgSpecFilt,3)
    imgSpecFiltSmooth(:,:,i)    =   wiener2(imgSpecFilt(:,:,i),filtSize);
end


%% Export final image for use in MESMA
imgSubFile  =   fullfile(dataDir,[siteDateCode 'SpecSub.tif']);
% Note that we have to reflip so that the exported images are in the correct
% orientation using the same R geo-referencing structure created when we
% imported the image in the first place
% geotiffwrite(imgSubFile,flipud(imgSpecSub),R,'CoordRefSysCode',epsgCode);
geotiffwrite(imgSubFile,flipud(imgSpecFiltSmooth),R,'CoordRefSysCode',epsgCode);
% Export a simple text file with wavelengths we can import to ENVI
writematrix(wvlSub,fullfile(dataDir,'wvlSubList.txt'));


%% Calculate Cyanobacterial index (CI) and display result
% Equation provided by Tyler King:
% CI=-Rrs(λ2)-Rrs(λ1)-(Rrs(λ3)-Rrs(λ1))*[(λ2-λ1)/(λ3-λ1)]
% λ1 = 662 nm λ2 = 681 nm λ3 = 709 nm
[~,iRedWvl1]    =   min(abs(wvlSub-662));
[~,iRedWvl2]    =   min(abs(wvlSub-681));
[~,iNirWvl]     =   min(abs(wvlSub-709));
ci      =   -1*(imgSpecFiltSmooth(:,:,iRedWvl2) - imgSpecFiltSmooth(:,:,iRedWvl1) - ...
            (imgSpecFiltSmooth(:,:,iNirWvl) - imgSpecFiltSmooth(:,:,iRedWvl1)) * ...
            (wvlSub(iRedWvl2)-wvlSub(iRedWvl1))/(wvlSub(iNirWvl)-wvlSub(iRedWvl1)));
ci      =   applyMask(ci,waterMask);
% Account for DESIS reflectance scale factor
ci      =   ci/10000;
figure
imagesc(ci,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
    'AlphaData',waterMask);
axis image; axis xy; axis on;
crameri('-bamako')
colorbar
title([siteDateCode ': Cyanobacterial Index'])
xlabel('Easting (m)')
ylabel('Northing (m)')


%% Export CI as Geo-TIFF
% Note that we have to reflip so that the exported image is in the correct
% orientation using the same R geo-referencing structure created when we
% imported the image in the first place
% Export CI image as a geotiff
maskFile=   fullfile(dataDir,[siteDateCode 'ci.tif']);
geotiffwrite(maskFile,flipud(ci),R,'CoordRefSysCode',epsgCode);

% 
% %% SS 665
% % Equation provided by Tyler King:
% % SS=-Rrs(λ2)-Rrs(λ1)-(Rrs(λ3)-Rrs(λ1))*[(λ2-λ1)/(λ3-λ1)]
% % λ1 = 620 nm λ2 = 665 nm λ3 = 681 nm
% [~,iRedWvl1]    =   min(abs(wvlSub-620));
% [~,iRedWvl2]    =   min(abs(wvlSub-665));
% [~,iNirWvl]     =   min(abs(wvlSub-681));
% % ss      =   -1*(imgSpecFiltSmooth(:,:,iRedWvl2) - imgSpecFiltSmooth(:,:,iRedWvl1) - ...
% %             (imgSpecFiltSmooth(:,:,iNirWvl) - imgSpecFiltSmooth(:,:,iRedWvl1)) * ...
% %             (wvlSub(iRedWvl2)-wvlSub(iRedWvl1))/(wvlSub(iNirWvl)-wvlSub(iRedWvl1)));
% ss      =   (imgSpecSub(:,:,iRedWvl2) - imgSpecSub(:,:,iRedWvl1) - ...
%             (imgSpecSub(:,:,iNirWvl) - imgSpecSub(:,:,iRedWvl1)) * ...
%             (wvlSub(iRedWvl2)-wvlSub(iRedWvl1))/(wvlSub(iNirWvl)-wvlSub(iRedWvl1)));
% ss      =   applyMask(ss,waterMask);
% % Account for DESIS reflectance scale factor
% ss      =   ss/10000;
% figure
% imagesc(ss,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
%     'AlphaData',waterMask);
% axis image; axis xy; axis on;
% % crameri('-bamako')
% colorbar
% % title([siteDateCode ': Cyanobacterial Index'])
% title('SS(665)')
% xlabel('Easting (m)')
% ylabel('Northing (m)')

