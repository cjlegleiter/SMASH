function [imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,...
            outputRoot,dataDir] = derSmash(smashMatFile,epsgCode) %#ok<STOUT>
% derSmash.m: Perform SMASH based on spectral derivatives, using outputs from
% initial run based on the original spectra
% 
%% derSmash.m:
%   Perform SMASH using derivative Spectra, based on an initial run using the
%   original spectra. This function repeats the multiple endmember spectral
%   mixture analysis using the output from an initial run based on the original
%   spectra but first applies the Savitzky-Golay smoothing filter and calculates
%   the special derivative for each pixel in an image before performing the new
%   MESMA run. The only user inputs required are the name of the *.mat file for
%   the original run and the updated EPSG code. The input parameters for the
%   smoothing filter are hardwired in the code. In addition, the full paths to
%   the derivative spectral library and derivative shade end number also are
%   hardwired. Outputs from this function are identical to those from
%   postSmash.m. In addition, a new set of figure files and a *.mat file with a
%   distinct name indicating the derivative are also generated.
% 
%% SYNTAX:
%   [imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,...
%             outputRoot,dataDir] = derSmash(smashMatFile,epsgCode);
% 
%% INPUTS:
%   smashMatFile:   Full file path to the *.mat file from the initial SMASH run 
%                   based on the original spectra. The variables within this
%                   file will be used within this function rather than
%                   recreating them with preSmash.m.
%   epsgCode:       EPSG code to be written to the geotiff file for the
%                   derivative image. This input is included to allow the code 
%                   to be updated to reflect the correct UTM zone.
% 
%% OUTPUTS:
%   outputRoot:     Root file name for MESMA outputs, with a time stamp used to
%                   identify the MESMA run     
%   imgDerFile:     String with full file path for the spectral derivative image
%   wvlDer:         Vector of wavelengths for the derivative image (one band
%                   less than the original)
%   R:              Spatial referencing structure array read from the *.mat file
%                   from the initial SMASH run based on the original image
%   siteCode:       Two-character site code read from the *.mat file from the 
%                   initial SMASH run based on the original image
%   dateCode:       Six-digit date code read from the *.mat file from the 
%                   initial SMASH run based on the original image
%   waterMask:      Binary water mask read from the *.mat file from the initial 
%                   SMASH run based on the original image
%   libDerFile:     String with full file path for derivative spectral library
% 
%% NOTES:
% > This function assumes that an initial run of SMASH based on the original
% Spectra has already been performed and can be used as input to this function.
% > Hardwired input parameters for the smoothing filter are a third order
% polynomial with a window size of seven and applying the filter twice.
%
%% FUNCTION SUMMARY:
%   [imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,...
%             outputRoot,dataDir] = derSmash(smashMatFile,epsgCode);
% 
%% CREDITS:
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 08/03/2021


%% MESMA OF DERIVATIVE SPECTRA
% See email correspondence on this topic, stimulated by UKL results
% Apply Savitzky-Golay smoothing filter to each spectrum in library and each
% pixel in image to reduce noise, then calculate spectral derivatives and apply
% MESMA to the first-derivative spectra


%% Use existing *.mat file from initial SMASH run using original spectra
% smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810.mat';
load(smashMatFile) %#ok<LOAD>
% Update EPSG code that we neglected to change from Lake Owasco
% epsgCode        =   32610; % For both Klamath and Detroit


%% Parameters for the Savitzky-Golay smoothing filter to be applied to raw data
nFilt   =   2;
order   =   3;
window  =   7;
% Compute spectral derivative as diff(R)/diff(lambda)
wvlDiff     =   diff(wvlSub');


%% Apply filter and calculate derivatives in a double for loop over the image pixels
imgDer      =   zeros(size(imgSpecSub,1),size(imgSpecSub,2),size(imgSpecSub,3)-1);
tmp         =   imgSpecSub;
disp('Computing spectral derivatives for image pixels ...')
for i = 1:size(tmp,1)
    for j = 1:size(tmp,2)
        tmpPix  =   squeeze(tmp(i,j,:));
        for iFilt=1:nFilt
            tmpPix  =   sgolayfilt(tmpPix,order,window);
        end    
        pixDiff         =   diff(tmpPix);
        pixDer          =   pixDiff./wvlDiff';
        imgDer(i,j,:)   =   pixDer;
    end
end


%% Write out the new derivative image as a geotiff
imgDerFile  =   fullfile(dataDir,[siteDateCode 'SpecDer.tif']);
% Note that we have to reflip so that the exported images are in the correct
% orientation using the same R geo-referencing structure created when we
% imported the image in the first place
geotiffwrite(imgDerFile,flipud(imgDer),R,'CoordRefSysCode',epsgCode);
% Export a simple text file with wavelengths we can import to ENVI
writematrix(wvlSub(1:end-1),fullfile(dataDir,'wvlDerList.txt')); %#ok<FNCOLND>


%% Now apply MESMA to derivative image using derivative library
libDerFile  =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibraryDerivative.sli';
shadeFile   =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterDerivative.sli';
disp('Performing MESMA based on derivative spectra ...')
outputRoot  =   runSMASH(dataDir,siteCode,dateCode,libDerFile,...
                          classNameField,imgDerFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
wvlDer      =   wvlSub(1:end-1);                       %#ok<FNCOLND>
                      
                       
%% Next call postSmash and save distinct figures and .mat file with derivative results
% close all
% disp('Post-processing MESMA output ...')
% [Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
%      taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
%                                        dateCode,waterMask,libDerFile,outputRoot);
% figDir      =   [dataDir '\figs'];
% if ~isfolder(figDir); mkdir(figDir); end                                   
% figure(1)
% set(gcf,'name',[siteDateCode '_RGBDerivative'])
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
% figure(2)
% set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
% figure(3)
% set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
% figure(4)
% set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
% figure(5)
% set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
% figure(6)
% set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
% figure(7)
% set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
% figure(8)
% set(gcf,'name',[siteDateCode '_RMSEDerivative'])
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
% saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
% 
% % Also save all results to a MAT file
% clear tmp imgSpecSub tmpPix ans i j
% save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))