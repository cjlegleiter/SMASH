%% Log for processing images from all sites via SMASH workflow
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 06/29/2021
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SMASHlog.m


%% Begin with the Lake Owasco data set we used as a prototype: OW20200823


%% OW20200823 pre-SMASH
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823';
siteCode    =   'OW';
dateCode    =   '20200823';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 18N is EPSG 32618;
epsgCode    =   32618;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823\rawImage\DESIS-HSI-L2A-DT0489826624_003-20200823T165204-V0210-SPECTRAL_IMAGE.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [-500 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
    

%% OW20200823 run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\AlgaeTaxaOnly.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% OW20200823 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))





%% On to another date for this site: OW20190805


%% OW20190805 pre-SMASH
close all
clc
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805';
siteCode    =   'OW';
dateCode    =   '20190805';
siteDateCode=   [siteCode dateCode];
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805\rawImage\DESIS-HSI-L2A-DT0348956444_003-20190805T185618-V0210-SPECTRAL_IMAGE.tif';
% Used [-500 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput);
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
    

%% OW20190805 run SMASH
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% OW20190805 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))





%% On to another date for this site: OW20190809


%% OW20190809 pre-SMASH
close all
clc
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190809';
siteCode    =   'OW';
dateCode    =   '20190809';
siteDateCode=   [siteCode dateCode];
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190809\rawImage\DESIS-HSI-L2A-DT0350460476_003-20190809T171514-V0210-SPECTRAL_IMAGE.tif';
% Used [-500 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput);
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
    

%% OW20190809 run SMASH
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% OW20190809 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))




%% On to another lake instead: Upper Klamath UK20200810

%% UK20200810 pre-SMASH
close all
clc
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810';
siteCode    =   'UK';
dateCode    =   '20200810';
siteDateCode=   [siteCode dateCode];
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\rawImage\DESIS-HSI-L2A-DT0485743928_009-20200810T185658-V0210-SPECTRAL_IMAGE.tif';
% Used [-350 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput);
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
    

%% UK20200810 run SMASH
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% UK20200810 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))





%% First date for Detroit - DL20190621 - THIS ONE DIDN'T WORK VERY WELL


%% DL20190621 pre-SMASH
close all
clc
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190621';
siteCode    =   'DL';
dateCode    =   '20190621';
siteDateCode=   [siteCode dateCode];
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190621\rawImage\DESIS-HSI-L2A-DT0332355592_002-20190621T212426-V0210-SPECTRAL_IMAGE_modified.tif';
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput);
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
    

%% DL20190621 run SMASH
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% DL20190621 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))



%% Next date for Detroit - DL20190819


%% DL20190819 pre-SMASH
close all
clc
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819';
siteCode    =   'DL';
dateCode    =   '20190819';
siteDateCode=   [siteCode dateCode];
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819\rawImage\DESIS-HSI-L2A-DT0354135644_002-20190819T215229-V0210-SPECTRAL_IMAGE.tif';
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput);
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
    

%% DL20190819 run SMASH
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% DL20190819 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))






%% Next date for Detroit - DL20200812


%% DL20200812 pre-SMASH
close all
clc
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200812';
siteCode    =   'DL';
dateCode    =   '20200812';
siteDateCode=   [siteCode dateCode];
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200812\rawImage\DESIS-HSI-L2A-DT0486486064_002-20200812T235111-V0210-SPECTRAL_IMAGE_modified.tif';
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput);
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
    

%% DL20200812 run SMASH
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% DL20200812 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))






%% Next date for Detroit - DL20200813


%% DL20200813 pre-SMASH
close all
clc
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813';
siteCode    =   'DL';
dateCode    =   '20200813';
siteDateCode=   [siteCode dateCode];
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813\rawImage\DESIS-HSI-L2A-DT0486825940_006-20200813T181126-V0210-SPECTRAL_IMAGE.tif';
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput);
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
    

%% DL20200813 run SMASH
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% DL20200813 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))






%% Next date for Detroit - DL20200828


%% DL20200828 pre-SMASH
close all
clc
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828';
siteCode    =   'DL';
dateCode    =   '20200828';
siteDateCode=   [siteCode dateCode];
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828\rawImage\DESIS-HSI-L2A-DT0492053448_002-20200828T173638-V0210-SPECTRAL_IMAGE.tif';
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput);
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
    

%% DL20200828 run SMASH
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% DL20200828 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))






%% NOTES FOR NEXT TIME - 7/1/2021
% Need to think about why some sites/dates produced no meaningful output - RMSE?
% Compare MESMA output to field samples in spreadsheet - which taxa are which?
% Circle back to simulated spectra and generate nEM X nEM matrix of fraction
% errors, similar to what we made for the NS3 score matrix; this will yield more
% insight on which taxa should be distinguishable




%% 7/1/2021: SIMULATED MESMA FRACTION ERRORS AS METRIC OF SEPARABILITY FOR ALL PAIRS OF TAXA
% Make use of the function we wrote for this purpose: HABsimSpec.m
% [outputRoot,fractions,inputMix,fracDiff,rmse,meanFracError,stdFracError,...
%           meanRmse,stdRmse] = HABsimSpec(libraryFile,wvl,iTaxa2mix1,iTaxa2mix2,...
%                                 iWater,mixCode,dataDir,minEmFrac,maxEmFrac,...
%                                 minShadeEmFrac,maxShadeEmFrac,maxRMSE)


%% Import spectral library and list of wavelengths, then plot the end-members
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\simSpec')
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\Algae+Water.sli';
% Get list of taxa from library and display numeric codes
Taxa    =   readtable([libraryFile(1:end-3) 'csv']);
Taxa    =   Taxa.Class;
TaxaNums=   cell(size(Taxa));
for i = 1:length(Taxa)
    TaxaNums{i} =   [num2str(i) ': ' Taxa{i}];
end
disp('Algal taxa numeric codes and taxa names:')
disp(TaxaNums)
hFig    =   uifigure('Position',[50 600 300 300],'Name','Algal Taxa Listing');
uilistbox(hFig,'Items',TaxaNums,'Position',[5 5 290 290]);
clear i h*

% Import the library
library =   envi2matlab(libraryFile);
% This is a 14 row by 158 column array with each end member as a row and each
% wavelength as a column. 
% Get list of wavelengths from a text file
wvlFile =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SimSpec\wvlSubList.txt';
wvl     =   load(wvlFile);


% Plot the spectra
figure
plot(wvl,library,'linewidth',2)
xlabel('Wavelength (nm)')
ylabel('Reflectance (%)')
title('Algal spectral library')
legend(TaxaNums)


%% Create arrays to store all pairwise summary metrics for taxa separability
SimSpec     =   struct('outputRoot',[],'fractions',[],'inputMix',[],'fracDiff',[],...
                   'rmse',[],'meanFracError',[],'stdFracError',[],'meanRmse',[],...
                   'stdRmse',[]);
% Don't count water as an end-member since it will be in all of the mixtures by
% definition and you can't include the same spectrum as two end-members or the
% results will be indeterminate and lead to NaN's
numEm       =   size(library,1)-1;
% Note that the fraction errors will need to be a 3D array because we will have
% a separate fraction error for each of the three end members that go into the
% mixture: the two algal taxa and water
meanFracErr =   zeros(numEm,numEm,3);
stdFracErr  =   zeros(numEm,numEm,3);
meanRmse    =   zeros(numEm,numEm);
stdRmse     =   zeros(numEm,numEm);
mixCode     =   cell(numEm,numEm);


%% Set up inputs we can recycle for all pairs
iWater  =   14;
dataDir =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SimSpec';
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


%% Loop over all taxa combinations, create simulated mixtures, perform MESMA, and summarize output
for i = 1:numEm
    for j = 1:numEm
        if i < j
            iTaxa2mix1  =   i;
            iTaxa2mix2  =   j;
            mixCode{i,j}=   [Taxa{i}(1:5) Taxa{j}(1:5)];
            disp(['Simulating spectral mixture for combination (' num2str(i) ',' ...
                  num2str(j) '): ' mixCode{i,j}])
            [SimSpec(i,j).outputRoot,SimSpec(i,j).fractions,...
             SimSpec(i,j).inputMix,SimSpec(i,j).fracDiff,SimSpec(i,j).rmse,...
             SimSpec(i,j).meanFracError,SimSpec(i,j).stdFracError,...
             SimSpec(i,j).meanRmse,SimSpec(i,j).stdRmse] ...
                        =   HABsimSpec(libraryFile,wvl,iTaxa2mix1,iTaxa2mix2,...
                                 iWater,mixCode{i,j},dataDir,minEmFrac,maxEmFrac,...
                                    minShadeEmFrac,maxShadeEmFrac,maxRMSE);
            meanFracErr(i,j,:)  =   SimSpec(i,j).meanFracError;
            stdFracErr(i,j,:)   =   SimSpec(i,j).stdFracError;
            meanRmse(i,j)       =   SimSpec(i,j).meanRmse;
            stdRmse(i,j)        =   SimSpec(i,j).stdRmse;
            close all
        end
    end
end


%% Save what we have so far
clear i j ans iTaxa* 
save HABsimSpecAllPairs.mat
% Delete all intermediate simulated mixture images and MESMA files in Windows
% explorer to save disk space


%% Visualize results: Mean fraction errors
% For now, average over the three end members fraction errors to get one value
% for each mixture
% Make a mask to only show the actual combinations
mask    =   zeros(numEm,numEm);
for i = 1:numEm
    for j = 1:numEm
        if i < j
            mask(i,j) = 1;
        end
    end
end
mask    =   logical(mask);
figure
imagesc(mean(meanFracErr,3),'alphadata',mask); 
axis equal; axis tight
crameri('imola'); colorbar
set(gca,'Xtick',1:numEm,'XtickLabel',Taxa(1:numEm))
set(gca,'Ytick',1:numEm,'YtickLabel',Taxa(1:numEm))
title('Mean MESMA fraction errors for algal taxa and water')
set(gcf,'name','SimSpecMeanFracErr')
saveas(gcf,get(gcf,'name'),'fig')
saveas(gcf,get(gcf,'name'),'jpg')


%% Visualize results: standard deviation of fraction errors
% For now, average over the three end members fraction errors to get one value
% for each mixture
figure
imagesc(mean(stdFracErr,3),'alphadata',mask); 
axis equal; axis tight
crameri('imola'); colorbar
set(gca,'Xtick',1:numEm,'XtickLabel',Taxa(1:numEm))
set(gca,'Ytick',1:numEm,'YtickLabel',Taxa(1:numEm))
title('Standard deviation of MESMA fraction errors for algal taxa and water')
set(gcf,'name','SimSpecStdFracErr')
saveas(gcf,get(gcf,'name'),'fig')
saveas(gcf,get(gcf,'name'),'jpg')


%% Visualize results: mean MESMA RMSE
figure
imagesc(meanRmse,'alphadata',mask); 
axis equal; axis tight
crameri('imola'); colorbar
set(gca,'Xtick',1:numEm,'XtickLabel',Taxa(1:numEm))
set(gca,'Ytick',1:numEm,'YtickLabel',Taxa(1:numEm))
title('Mean MESMA RMSE for algal taxa and water')
set(gcf,'name','SimSpecMeanRMSE')
saveas(gcf,get(gcf,'name'),'fig')
saveas(gcf,get(gcf,'name'),'jpg')


%% Visualize results: standard deviation of MESMA RMSE
figure
imagesc(stdRmse,'alphadata',mask); 
axis equal; axis tight
crameri('imola'); colorbar
set(gca,'Xtick',1:numEm,'XtickLabel',Taxa(1:numEm))
set(gca,'Ytick',1:numEm,'YtickLabel',Taxa(1:numEm))
title('Standard deviation of MESMA RMSE for algal taxa and water')
set(gcf,'name','SimSpecStdRMSE')
saveas(gcf,get(gcf,'name'),'fig')
saveas(gcf,get(gcf,'name'),'jpg')







%% 7/13/2021: REPEAT UKL EXCLUDING TOLYPOTHRIX FROM LIBRARY
% Use the subset library we just created:
%   C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySubset.sli
% Work on making a general function for subsetting the library


%% ALSO NEXT TIME: MESMA OF DERIVATIVE SPECTRA
% See email correspondence on this topic, stimulated by UKL results
% Apply Savitzky-Golay smoothing filter to each spectrum in library and each
% pixel in image to reduce noise, then calculate spectral derivatives and apply
% MESMA to the first-derivative spectra




%% **** 7/15/2021: REPEAT UKL EXCLUDING TOLYPOTHRIX FROM LIBRARY
% Use the new function subsetLibrary.m
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\Algae+Water.sli';
librarySub      =   subsetLibrary(libraryFile);
% Select taxa to include: [1:12]
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\LibrarySubset.sli';


%% UK20200810 pre-SMASH
close all
clc
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810';
siteCode    =   'UK';
dateCode    =   '20200810';
siteDateCode=   [siteCode dateCode];
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\rawImage\DESIS-HSI-L2A-DT0485743928_009-20200810T185658-V0210-SPECTRAL_IMAGE.tif';
epsgCode    =   32610;
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [-350 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput);
imtool close all
close(1)


%% UK20200810 run SMASH
% The rest should be consistent for all sites
% Metadata field with end-member names
classNameField  =   'Class';
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

outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% UK20200810 post-SMASH - save distinct figures and a .mat file with new no-Toly results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot); 
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBNoToly'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationNoToly'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramNoToly'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapNoToly'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsNoToly'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsNoToly'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBNoToly'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSENoToly'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode 'NoToly.mat']))


%% Looks like this worked OK
% We have no Toly when we exclude this from the library and more AFA. Added the
% figures to the PowerPoint presentation for smash.




%% **** 7/29/2021: Adapt smash for spectral derivatives
% MESMA OF DERIVATIVE SPECTRA
% See email correspondence on this topic, stimulated by UKL results
% Apply Savitzky-Golay smoothing filter to each spectrum in library and each
% pixel in image to reduce noise, then calculate spectral derivatives and apply
% MESMA to the first-derivative spectra


%% Apply Savitzky-Golay smoothing filter to each spectrum in library
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\Algae+Water.sli';
Taxa    =   readtable([libraryFile(1:end-3) 'csv']);
Taxa    =   Taxa.Class;
TaxaNums=   cell(size(Taxa));
for i = 1:length(Taxa)
    TaxaNums{i} =   [num2str(i) ': ' Taxa{i}];
end
disp('Algal taxa numeric codes and taxa names:')
disp(TaxaNums)
hFig    =   uifigure('Position',[50 600 300 300],'Name','Algal Taxa Listing');
hLabel  =   uilistbox(hFig,'Items',TaxaNums,'Position',[5 5 290 290]);

% Actual spectral library data file hard-wired into code
% We had to make the spectral library into an image data file to bring it into
% MATLAB
mergedLibraryFile   =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\Algae+WaterLibraryAsImage.dat';
library             =   envi2matlab(mergedLibraryFile);
% This is a 14 row by 158 column array with each end member as a row and each
% wavelength as a column. 
% Bring in list of wavelengths from our OWASCO .mat file
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SMASHowasco.mat','wvlSub');
% Plot the spectra
figure
plot(wvlSub,library)
xlabel('Wavelength (nm)')
ylabel('Reflectance (%)')
title('Original algae and water spectral library')
legend(TaxaNums)

% Apply smoothing filter to spectra
% Parameters for the Savitzky-Golay smoothing filter to be applied to raw data:
nFilt   =   2;
order   =   3;
window  =   7;
tmp     =   library';   
for iFilt=1:nFilt
    tmp =   sgolayfilt(tmp,order,window);
end    
libFilt =   tmp';
figure
plot(wvlSub,libFilt)
xlabel('Wavelength (nm)')
ylabel('Reflectance (%)')
title('Smoothed algae and water spectral library')
legend(TaxaNums)


%% Compute spectral derivative as diff(R)/diff(lambda)
libDiff     =   diff(libFilt,1,2);
wvlDiff     =   diff(wvlSub');
libDer      =   libDiff./wvlDiff;
figure
plot(wvlSub(1:end-1),libDer)
xlabel('Wavelength (nm)')
ylabel('Reflectance (%)')
title('First derivative algae and water spectral library')
legend(TaxaNums)


%% Write put the derivative spectra as a new library file
matlab2envi(libDer(1:13,:),'LibraryDerivative.sli',[],[],[],wvlSub(1:end-1));
Header      =   specLibHeader(libDer(1:13,:),'LibraryDerivative.sli',Taxa(1:13),wvlSub(1:end-1)); %#ok<*NASGU>
% Also make the corresponding csv metadata file
Metadata    =   table(Taxa(1:13),Taxa(1:13),'VariableNames',{'Name','Class'});
writetable(Metadata,'LibraryDerivative.csv')


%% Also make a new derivative shade file
matlab2envi(libDer(14,:),'WaterDerivative.sli',[],[],[],wvlSub(1:end-1));
Header      =   specLibHeader(libDer(14,:),'WaterDerivative.sli',Taxa(14),wvlSub(1:end-1)); %#ok<*NASGU>
% Also make the corresponding csv metadata file
Metadata    =   table(Taxa(14),Taxa(14),'VariableNames',{'Name','Class'});
writetable(Metadata,'WaterDerivative.csv')


%% Now apply this approach to the image in a double for loop over the pixels
imgDer      =   zeros(size(imgSpecSub,1),size(imgSpecSub,2),size(imgSpecSub,3)-1);
tmp         =   imgSpecSub;
for i = 1:size(tmp,1)
    for j = 1:size(tmp,2)
        disp(['Computing derivative spectrum for pixel in row ' num2str(i) ...
              ', column ' num2str(j) ' ...'])
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
writematrix(wvlSub(1:end-1),fullfile(dataDir,'wvlDerList.txt'));


%% Now apply MESMA to derivative image using derivative library
libDerFile  =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\LibraryDerivative.sli';
shadeFile   =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\WaterDerivative.sli';
outputRoot  =   runSMASH(dataDir,siteCode,dateCode,libDerFile,...
                          classNameField,imgDerFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% UK20200810 post-SMASH - save distinct figures and a .mat file with new derivative results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlSub(1:end-1),R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot); 
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))


%% Looks like this worked OK
% We have no Toly when we use spectral derivatives as input to the MESMA and AFA
% is the dominant taxon, consistent with the field samples. Added the figures to
% the PowerPoint presentation for SMASH.










%% **** 8/3/2021: GENERALIZE SPECTRAL DERIVATIVE CODE AND APPLY TO OTHER LAKES
% See the new function derSmash.m, which we can base on the UKL example above


%% Test derSmash.m using UKL inputs
% Use existing *.mat file from initial SMASH run using original spectra
smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810.mat';
% Update EPSG code that we neglected to change from Lake Owasco
epsgCode        =   32610; % For both Klamath and Detroit
[imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,outputRoot,dataDir] ...
                =   derSmash(smashMatFile,epsgCode);
            
            
%% Next call postSmash and save distinct figures and .mat file with derivative results
close all
disp('Post-processing MESMA output ...')
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot);
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear tmp imgSpecSub tmpPix ans i j
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))


%% Looks like this will work, so now apply this template to other sites/dates ...


%% OW20190805
% Use existing *.mat file from initial SMASH run using original spectra
smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805\OW20190805.mat';
% Update EPSG code that we neglected to change from Lake Owasco
epsgCode        =   32618; % For Owasco
[imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,outputRoot,dataDir] ...
                =   derSmash(smashMatFile,epsgCode);

            
%% Post-SMASH            
close all
disp('Post-processing MESMA output ...')
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot);
figDir      =   [dataDir '\figs'];
siteDateCode=   [siteCode dateCode];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear tmp imgSpecSub tmpPix ans i j
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))


%% OW20190809
% Use existing *.mat file from initial SMASH run using original spectra
smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190809\OW20190809.mat';
% Update EPSG code that we neglected to change from Lake Owasco
epsgCode        =   32618; % For Owasco
[imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,outputRoot,dataDir] ...
                =   derSmash(smashMatFile,epsgCode);

            
%% Post-SMASH            
close all
disp('Post-processing MESMA output ...')
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot);
figDir      =   [dataDir '\figs'];
siteDateCode=   [siteCode dateCode];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear tmp imgSpecSub tmpPix ans i j
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))


%% OW20200823
% Use existing *.mat file from initial SMASH run using original spectra
smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823\OW20200823.mat';
% Update EPSG code that we neglected to change from Lake Owasco
epsgCode        =   32618; % For Owasco
[imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,outputRoot,dataDir] ...
                =   derSmash(smashMatFile,epsgCode);

            
%% Post-SMASH            
close all
disp('Post-processing MESMA output ...')
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot);
figDir      =   [dataDir '\figs'];
siteDateCode=   [siteCode dateCode];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear tmp imgSpecSub tmpPix ans i j
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))


%% DL20190621
% Use existing *.mat file from initial SMASH run using original spectra
smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190621\DL20190621.mat';
% Update EPSG code that we neglected to change from Lake Owasco
epsgCode        =   32618; % For Owasco
[imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,outputRoot,dataDir] ...
                =   derSmash(smashMatFile,epsgCode);

            
%% Post-SMASH            
close all
disp('Post-processing MESMA output ...')
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot);
figDir      =   [dataDir '\figs'];
siteDateCode=   [siteCode dateCode];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear tmp imgSpecSub tmpPix ans i j
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))


%% DL20190819
% Use existing *.mat file from initial SMASH run using original spectra
smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819\DL20190819.mat';
% Update EPSG code that we neglected to change from Lake Owasco
epsgCode        =   32618; % For Owasco
[imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,outputRoot,dataDir] ...
                =   derSmash(smashMatFile,epsgCode);

            
%% Post-SMASH            
close all
disp('Post-processing MESMA output ...')
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot);
figDir      =   [dataDir '\figs'];
siteDateCode=   [siteCode dateCode];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear tmp imgSpecSub tmpPix ans i j
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))


%% DL20200812
% Use existing *.mat file from initial SMASH run using original spectra
smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200812\DL20200812.mat';
% Update EPSG code that we neglected to change from Lake Owasco
epsgCode        =   32618; % For Owasco
[imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,outputRoot,dataDir] ...
                =   derSmash(smashMatFile,epsgCode);

            
%% Post-SMASH            
close all
disp('Post-processing MESMA output ...')
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot);
figDir      =   [dataDir '\figs'];
siteDateCode=   [siteCode dateCode];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear tmp imgSpecSub tmpPix ans i j
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))


%% DL20200813
% Use existing *.mat file from initial SMASH run using original spectra
smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813\DL20200813.mat';
% Update EPSG code that we neglected to change from Lake Owasco
epsgCode        =   32610; % For Detroit
[imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,outputRoot,dataDir] ...
                =   derSmash(smashMatFile,epsgCode);

            
%% Post-SMASH            
close all
disp('Post-processing MESMA output ...')
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot);
figDir      =   [dataDir '\figs'];
siteDateCode=   [siteCode dateCode];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear tmp imgSpecSub tmpPix ans i j
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))


%% DL20200828
% Use existing *.mat file from initial SMASH run using original spectra
smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828\DL20200828.mat';
% Update EPSG code that we neglected to change from Lake Owasco
epsgCode        =   32610; % For Detroit
[imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,outputRoot,dataDir] ...
                =   derSmash(smashMatFile,epsgCode);

            
%% Post-SMASH            
close all
disp('Post-processing MESMA output ...')
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot);
figDir      =   [dataDir '\figs'];
siteDateCode=   [siteCode dateCode];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear tmp imgSpecSub tmpPix ans i j
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))












%% **** 8/4/2021: PROCESS NEW DATA FROM GRAPEVINE LAKE, TEXAS
% Email from Kurt Carpenter on field samples
% The dominant taxa in Grapevine on 8/19/2020 were:
% 
% Raphidiopsis, aka Cylindrospermopsis (note this is distinct from
% Cylindrospermum in our library).The MUM is unique! (Right Adam!?) 🙂 Other
% pseudonyms include Cylindrospermopsis philippinensis, Anabaenopsis
% raciborskii, and Cylindrospermopsis raciborskii
% 
% Pseudanabaena coming in second is interesting, this cyano is often found along
% with Microcystis, clearly visible in the MC mucilage (see pic below) they are
% the short linear filaments among a sea of round Microcystis cells. So we might
% also see MC in the imagery.
% 
% Cyclostephanos is a diatom Peridinium is a dinoflagellate; some of these
% organisms produce toxins, but low percent of total biovolume
% 
% 
% 	relative_total_biovolume_percent
% Raphidiopsis	24%
% Pseudanabaena	20%
% Cyclostephanos	8%
% Peridinium	7%


%% GV20200813 pre-SMASH
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813';
siteCode    =   'GV';
dateCode    =   '20200813';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 14S is EPSG 32614;
epsgCode    =   32614;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813\rawImage\DESIS-HSI-L2A-DT0486845936_004-20200813T230925-V0213-SPECTRAL_IMAGE.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [-500 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
    

%% GV20200813 run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\AlgaeTaxaOnly.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% GV20200813 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))


%% Second tile for this site


%% GV20200813b pre-SMASH
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813b';
siteCode    =   'GV';
dateCode    =   '20200813b';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 14S is EPSG 32614;
epsgCode    =   32614;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813b\rawImage\DESIS-HSI-L2A-DT0486845936_005-20200813T230925-V0213-SPECTRAL_IMAGE.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [-500 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
    

%% GV20200813b run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\AlgaeTaxaOnly.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% GV20200813b post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))


%% Also apply to spectral derivative images for these two tiles ...


%% GV20200813
% Use existing *.mat file from initial SMASH run using original spectra
smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813\GV20200813.mat';
% Update EPSG code that we neglected to change from Lake Owasco
epsgCode        =   32614;
[imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,outputRoot,dataDir] ...
                =   derSmash(smashMatFile,epsgCode);

            
%% Post-SMASH            
close all
disp('Post-processing MESMA output ...')
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot);
figDir      =   [dataDir '\figs'];
siteDateCode=   [siteCode dateCode];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear tmp imgSpecSub tmpPix ans i j
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))


%% GV20200813b
% Use existing *.mat file from initial SMASH run using original spectra
smashMatFile    =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813b\GV20200813b.mat';
% Update EPSG code that we neglected to change from Lake Owasco
epsgCode        =   32614;
[imgDerFile,wvlDer,R,siteCode,dateCode,waterMask,libDerFile,outputRoot,dataDir] ...
                =   derSmash(smashMatFile,epsgCode);

            
%% Post-SMASH            
close all
disp('Post-processing MESMA output ...')
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgDerFile,wvlDer,R,siteCode,...
                                       dateCode,waterMask,libDerFile,outputRoot);
figDir      =   [dataDir '\figs'];
siteDateCode=   [siteCode dateCode];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEDerivative'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear tmp imgSpecSub tmpPix ans i j
save(fullfile(dataDir,[siteDateCode 'Derivative.mat']))









%% **** 8/10/2021: REPEAT UKL WITH SPECTRAL SMOOTHING BUT NOT TAKING DERIVATIVE, AND NO TOLY


%% Load *.mat file with results for UKL from initial SMASH run
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810.mat')


%% Apply Savitzky-Golay smoothing filter to each spectrum in subset library
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySubset.sli';
Taxa    =   readtable([libraryFile(1:end-3) 'csv']);
Taxa    =   Taxa.Class;
TaxaNums=   cell(size(Taxa));
for i = 1:length(Taxa)
    TaxaNums{i} =   [num2str(i) ': ' Taxa{i}];
end
disp('Algal taxa numeric codes and taxa names:')
disp(TaxaNums)
hFig    =   uifigure('Position',[50 600 300 300],'Name','Algal Taxa Listing');
hLabel  =   uilistbox(hFig,'Items',TaxaNums,'Position',[5 5 290 290]);

% Actual spectral library data file hard-wired into code
% We had to make the spectral library into an image data file to bring it into
% MATLAB
mergedLibraryFile   =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\Algae+WaterLibraryAsImage.dat';
library             =   envi2matlab(mergedLibraryFile);
% Remove Toly (row 13)
library(13,:)       =   [];
% This is now a 13 row by 158 column array with each end member as a row and each
% wavelength as a column. 
% Bring in list of wavelengths from our OWASCO .mat file
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SMASHowasco.mat','wvlSub');
% Plot the spectra
figure
plot(wvlSub,library(1:12,:))
xlabel('Wavelength (nm)')
ylabel('Reflectance (%)')
title('Original algae spectral library')
legend(TaxaNums)

% Apply smoothing filter to spectra
% Parameters for the Savitzky-Golay smoothing filter to be applied to raw data:
nFilt   =   2;
order   =   3;
window  =   7;
tmp     =   library';   
for iFilt=1:nFilt
    tmp =   sgolayfilt(tmp,order,window);
end    
libFilt =   tmp';
figure
plot(wvlSub,libFilt(1:12,:))
xlabel('Wavelength (nm)')
ylabel('Reflectance (%)')
title('Smoothed algae spectral library')
legend(TaxaNums)


%% Write out the smoothed spectra as a new library file
matlab2envi(libFilt(1:12,:),'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli',...
                [],[],[],wvlSub);
Header      =   specLibHeader(libFilt(1:12,:),...
                'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli',...
                Taxa(1:12),wvlSub); %#ok<*NASGU>
% Also make the corresponding csv metadata file
Metadata    =   table(Taxa(1:12),Taxa(1:12),'VariableNames',{'Name','Class'});
writetable(Metadata,'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.csv')


%% Also make a new smoothed shade file
matlab2envi(libFilt(13,:),...
    'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli',...
    [],[],[],wvlSub);
Header      =   specLibHeader(libFilt(13,:),...
                'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli',...
               {'Water'},wvlSub); %#ok<*NASGU>
% Also make the corresponding csv metadata file
Metadata    =   table({'Water'},{'Water'},'VariableNames',{'Name','Class'});
writetable(Metadata,'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.csv')


%% Now apply this approach to the image in a double for loop over the pixels
imgSpecFilt =   zeros(size(imgSpecSub,1),size(imgSpecSub,2),size(imgSpecSub,3));
tmp         =   imgSpecSub;
for i = 1:size(tmp,1)
    for j = 1:size(tmp,2)
        disp(['Smoothing spectrum for pixel in row ' num2str(i) ...
              ', column ' num2str(j) ' ...'])
        tmpPix  =   squeeze(tmp(i,j,:));
        for iFilt=1:nFilt
            tmpPix  =   sgolayfilt(tmpPix,order,window);
        end    
        imgSpecFilt(i,j,:)   =   tmpPix;
    end
end


%% Write out the new spectrally smoothed image as a geotiff
imgSpecFiltFile  =   fullfile(dataDir,[siteDateCode 'SpecFilt.tif']);
% Note that we have to reflip so that the exported images are in the correct
% orientation using the same R geo-referencing structure created when we
% imported the image in the first place
epsgCode        =   32610;
geotiffwrite(imgSpecFiltFile,flipud(imgSpecFilt),R,'CoordRefSysCode',epsgCode);


%% Now apply MESMA to spectrally smoothed image using smoothed library
libFiltFile =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
shadeFile   =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';
outputRoot  =   runSMASH(dataDir,siteCode,dateCode,libFiltFile,...
                          classNameField,imgSpecFiltFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% UK20200810 post-SMASH - save distinct figures and a .mat file with new spectrally smoothed results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSpecFiltFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libFiltFile,outputRoot); 
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBSmoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationSmoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramSmoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapSmoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsSmoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsSmoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBSmoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSESmoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode 'Smoothed.mat']))


%% Looks like this worked OK, but not a big difference relative to initial run
% Main chnge relative to the run with no Toly but the original spectra was a
% less pixelated and noisy spatial pattern, which is nice, at least cosmetically



%% Try applying a spatial smoothing filter as well
% doc wiener2
filtSize            =   [3 3];
imgSpecFiltSmooth   =   zeros(size(imgSpecFilt));
for i = 1:size(imgSpecFilt,3)
    imgSpecFiltSmooth(:,:,i)    =   wiener2(imgSpecFilt(:,:,i),filtSize);
end
% Compare images for a band at 675.8 nm
iDispWvl            =   108;
figure
tiledlayout(1,2)
hAxOriginal =   nexttile;
% Set up 98% contrast stretch, excluding the zeros
tmp         =   imgSpecFilt(:,:,iDispWvl);
bad         =   isnan(tmp) | tmp==0;
tmp(bad)    =   nan;
dispRange   =   prctile(tmp(:),[2.5 97.5]);
imshow(imgSpecFilt(:,:,iDispWvl),...
       'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
       'displayrange',dispRange);
colormap gray; axis image; axis xy; axis on
title('Spectrally filtered Image')
xlabel('Easting (m)')
ylabel('Northing (m)')
hAxFilt     =   nexttile;
imshow(imgSpecFiltSmooth(:,:,iDispWvl),...
       'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
       'displayrange',dispRange);
colormap gray; axis image; axis xy; axis on
title('Spectrally Filtered and Spatially Smoothed Image')
xlabel('Easting (m)')
ylabel('Northing (m)')
linkaxes([hAxOriginal hAxFilt])


%% Write out the new spectrally filtered and spatially smoothed image as a geotiff
imgSpecFiltSmoothFile  =   fullfile(dataDir,[siteDateCode 'SpecFiltSmooth.tif']);
% Note that we have to reflip so that the exported images are in the correct
% orientation using the same R geo-referencing structure created when we
% imported the image in the first place
epsgCode        =   32610;
geotiffwrite(imgSpecFiltSmoothFile,flipud(imgSpecFiltSmooth),R,'CoordRefSysCode',epsgCode);


%% Now apply MESMA to spectrally filtered and spatially smoothed image using filtered library
libFiltFile =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
shadeFile   =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';
outputRoot  =   runSMASH(dataDir,siteCode,dateCode,libFiltFile,...
                          classNameField,imgSpecFiltSmoothFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% UK20200810 post-SMASH - save distinct figures and a .mat file with new spectrally filtered and spatially smoothed results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSpecFiltSmoothFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libFiltFile,outputRoot); 
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGBFiltered+Smoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassificationFiltered+Smoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogramFiltered+Smoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMapFiltered+Smoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractionsFiltered+Smoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractionsFiltered+Smoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGBFiltered+Smoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSEFiltered+Smoothed'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
clear imgSpecFilt
save(fullfile(dataDir,[siteDateCode 'Filtered+Smoothed.mat']))


%% Looks like this worked OK, but not a big difference relative to initial run
% Main chnge relative to the run with no Toly but the original spectra was a
% less pixelated and noisy spatial pattern, which is nice, at least cosmetically




%% **** 8/10/2021: ESTIMATE MEAN RMSE FOR PIXELS CLASSIFIED AS EACH TAXA
% % Get list of unique taxa ID numbers from classified image
% uniTaxa     =   unique(modelSum(~isnan(modelSum))); % Just 1-12 in this case
meanRmse    =   zeros(size(Taxa));
for i = 1:length(Taxa)
   inTaxa       =   find(modelSum == i & waterMask);
   tmpRmse      =   rmse(inTaxa);
   meanRmse(i)  =   mean(tmpRmse,'omitnan');
end
figure
bar(meanRmse);
ylabel('RMSE')
set(gca,'xticklabels',Taxa)
title([siteDateCode ': MESMA RMSE by algal taxon'])






%% **** 8/10/2021: EXTRACT FRACTIONS FOR A GIVEN POINT SAMPLE
% Rattlesnake Point in UKL: 42.344214	-121.858296
fieldLat        =   42.344214;
fieldLon        =   -121.858296;
% Convert to UTM 
[fieldX,fieldY] =   projfwd(projcrs(epsgCode),fieldLat,fieldLon);
% Extract fractions for all taxa at this location
fracAtField     =   zeros(length(Taxa)+1,1);
for i = 1:length(fracAtField)
    tmp             =   impixel(R.XWorldLimits,R.YWorldLimits,...
                                fractions(:,:,i),fieldX,fieldY);
    fracAtField(i)  =   tmp(1);
end
figure
bar(fracAtField)
ylabel('End member fraction')
set(gca,'xticklabels',[Taxa; {'Water'}])
title([siteDateCode ': end member fractions for algal taxa'])


%% **** 8/18/2021: CALCULATE CYANOBACTERIAL INDEX AS PART OF preSmash.m
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
figure
imagesc(ci,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
    'AlphaData',waterMask);
axis image; axis xy; axis on;
crameri('-bamako')
colorbar
title([siteDateCode ': Cyanobacterial Index'])
xlabel('Easting (m)')
ylabel('Northing (m)')










%%
%%
%% **** 9/24/2021: REPROCESS ALL IMAGES WITH COMMON WORKFLOW AND LIBRARY


%% GV20200813 pre-SMASH
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813';
siteCode    =   'GV';
dateCode    =   '20200813';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 14S is EPSG 32614;
epsgCode    =   32614;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813\rawImage\DESIS-HSI-L2A-DT0486845936_004-20200813T230925-V0213-SPECTRAL_IMAGE.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [0 1500] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,ci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
figure(2)
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_CI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');


%% GV20200813 run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% GV20200813 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))


%% Second tile for this site


%% GV20200813b pre-SMASH
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813b';
siteCode    =   'GV';
dateCode    =   '20200813b';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 14S is EPSG 32614;
epsgCode    =   32614;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813b\rawImage\DESIS-HSI-L2A-DT0486845936_005-20200813T230925-V0213-SPECTRAL_IMAGE.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [0 1500] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,ci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
figure(2)
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_CI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');


%% GV20200813b run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% GV20200813b post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))






%% OW20190805 pre-SMASH
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805';
siteCode    =   'OW';
dateCode    =   '20190805';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 18T is EPSG 32618;
epsgCode    =   32618;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805\rawImage\DESIS-HSI-L2A-DT0348956444_003-20190805T185618-V0210-SPECTRAL_IMAGE.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [0 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,ci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
figure(2)
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_CI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');


%% OW20190805 run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% OW20190805 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))







%% OW20200823 pre-SMASH
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823';
siteCode    =   'OW';
dateCode    =   '20200823';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 18T is EPSG 32618;
epsgCode    =   32618;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823\rawImage\DESIS-HSI-L2A-DT0489826624_003-20200823T165204-V0210-SPECTRAL_IMAGE.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [0 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,ci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
figure(2)
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_CI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');


%% OW20200823 run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% OW20200823 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))







%% DL20190819 pre-SMASH
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819';
siteCode    =   'DL';
dateCode    =   '20190819';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 10T is EPSG 32610;
epsgCode    =   32610;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819\rawImage\DESIS-HSI-L2A-DT0354135644_002-20190819T215229-V0210-SPECTRAL_IMAGE.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [0 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,ci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
figure(2)
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_CI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');


%% DL20190819 run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% DL20190819 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))







%% DL20200813 pre-SMASH
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813';
siteCode    =   'DL';
dateCode    =   '20200813';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 10T is EPSG 32610;
epsgCode    =   32610;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813\rawImage\DESIS-HSI-L2A-DT0486825940_006-20200813T181126-V0210-SPECTRAL_IMAGE.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [0 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,ci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
figure(2)
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_CI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');


%% DL20200813 run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% DL20200813 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))







%% DL20200828 pre-SMASH
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828';
siteCode    =   'DL';
dateCode    =   '20200828';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 10T is EPSG 32610;
epsgCode    =   32610;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828\rawImage\DESIS-HSI-L2A-DT0492053448_002-20200828T173638-V0210-SPECTRAL_IMAGE.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [0 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,ci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
figure(2)
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_CI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');


%% DL20200828 run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% DL20200828 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))







%% UK20200810 pre-SMASH
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810';
siteCode    =   'UK';
dateCode    =   '20200828';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 10T is EPSG 32610;
epsgCode    =   32610;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\rawImage\DESIS-HSI-L2A-DT0485743928_009-20200810T185658-V0210-SPECTRAL_IMAGE.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [-350 1000] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,ci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
figure(2)
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_CI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');


%% UK20200810 run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
                      
%% UK20200810 post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))






%%
%%
%% **** 10/4/2021: REPEAT SPECTRAL SEPARABILITY ANALYSIS WITH SMOOTHED LIBRARY SUBSET
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')

%%
%% Import smoothed algal library and water (shade) end member
% The MATLAB Hyperspectral Imaging Toolbox has several spectra for this purpose
%   web(fullfile(docroot, 'images/hyperspectral-image-processing.html?s_tid=CRUX_lftnav'))
% Import the library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
waterFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';
% Get list of taxa from library and display numeric codes
Taxa    =   readtable([libraryFile(1:end-3) 'csv']);
Taxa    =   Taxa.Class;
Taxa    =   [Taxa; {'Water'}];
TaxaNums=   cell(size(Taxa));
for i = 1:length(Taxa)
    TaxaNums{i} =   [num2str(i) ': ' Taxa{i}];
end
disp('Algal taxa numeric codes and taxa names:')
disp(TaxaNums)
hFig    =   uifigure('Position',[50 600 300 300],'Name','Algal Taxa Listing');
hLabel  =   uilistbox(hFig,'Items',TaxaNums,'Position',[5 5 290 290]);
% We had to make the spectral library into an image data file to bring it into
% MATLAB
library =   envi2matlab(libraryFile);
% This is a 12 row by 158 column array with each end member as a row and each
% wavelength as a column. 
water   =   envi2matlab(waterFile);
% This one just has water, so merge
libWater=   [library; water];
% Bring in list of wavelengths from our OWASCO .mat file
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SMASHowasco.mat','wvlSub');

%%
%% Plot the spectra
figure
h       =   plot(wvlSub,libWater,'linewidth',2);
TaxaLegend  =   Taxa;
for i = 1:length(Taxa)-1
    TaxaLegend{i}   =   ['{\it{' TaxaLegend{i} '}}'];
end
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
xlabel('Wavelength (nm)')
ylabel('Reflectance (%)')
title('Smoothed, merged algae and water spectral library')
% legend(TaxaNums)
legend(TaxaLegend)
set(gcf,'name','SmoothedLibrySpectra+Water')
saveas(gcf,fullfile('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library',get(gcf,'name')),'fig')

%%
%% Metrics of separability in a numEm X numEm matrix
% The Hyperspectral Imaging Toolbox has several measures of spectral similarity
% available but the best summary might be the normalized spectral similiarity
% score (NS3), which combines both the overall brightness (Euclidean distance)
% and spectral shape (spectral angle) in a summary score. See documentation for
% the function ns3: "The test spectrum with the minimum score value indicates
% highest similarity to the reference endmember. On the other hand, the test
% spectrum with the maximum score value has the highest spectral variability."
numEm   =   size(libWater,1);
ns3score=   zeros(numEm,numEm);
for i = 1:numEm
    for j = 1:numEm
        if i < j
            ns3score(i,j)   =   ns3(libWater(i,:),libWater(j,:));
        end
    end
end
figure
imagesc(ns3score,'alphadata',logical(ns3score)); 
axis equal; axis tight
crameri('imola'); colorbar
% set(gca,'Xtick',1:14,'XtickLabel',Taxa)
% set(gca,'Ytick',1:14,'YtickLabel',Taxa)
set(gca,'Xtick',1:14,'XtickLabel',TaxaLegend)
set(gca,'Ytick',1:14,'YtickLabel',TaxaLegend)
title('Normalized spectral separability scores for algal taxa and water')
% Looks like anabena is relatively distinct from most of the other taxa by this
% metric, so let's try making a mixture of anabena vs. microcystis, which is
% actually harmful and has been observed in Lake Owasco
clear i j ans h*
set(gcf,'name','SpectralSeparabilityScoreMatrixPlot')
saveas(gcf,fullfile('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\',get(gcf,'name')),'fig')

%%
%% 10/4/2021: SIMULATED MESMA FRACTION ERRORS AS METRIC OF SEPARABILITY FOR ALL PAIRS OF TAXA
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SimSpec')
% Make use of the function we wrote for this purpose: HABsimSpec.m
% [outputRoot,fractions,inputMix,fracDiff,rmse,meanFracError,stdFracError,...
%           meanRmse,stdRmse] = HABsimSpec(libraryFile,wvl,iTaxa2mix1,iTaxa2mix2,...
%                                 iWater,mixCode,dataDir,minEmFrac,maxEmFrac,...
%                                 minShadeEmFrac,maxShadeEmFrac,maxRMSE)

%%
%% Create arrays to store all pairwise summary metrics for taxa separability
SimSpec     =   struct('outputRoot',[],'fractions',[],'inputMix',[],'fracDiff',[],...
                   'rmse',[],'meanFracError',[],'stdFracError',[],'meanRmse',[],...
                   'stdRmse',[]);
% Don't count water as an end-member since it will be in all of the mixtures by
% definition and you can't include the same spectrum as two end-members or the
% results will be indeterminate and lead to NaN's
numEm       =   size(libWater,1)-1;
% Note that the fraction errors will need to be a 3D array because we will have
% a separate fraction error for each of the three end members that go into the
% mixture: the two algal taxa and water
meanFracErr =   zeros(numEm,numEm,3);
stdFracErr  =   zeros(numEm,numEm,3);
meanRmse    =   zeros(numEm,numEm);
stdRmse     =   zeros(numEm,numEm);
mixCode     =   cell(numEm,numEm);

%%
%% Set up inputs we can recycle for all pairs
iWater  =   13;
dataDir =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SimSpec';
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

%%
%% Anabena-Microcystis test case
i           =   1;
j           =   9;
wvl         =   wvlSub;     
iTaxa2mix1  =   i;
iTaxa2mix2  =   j;
mixCode{i,j}=   [Taxa{i}(1:5) Taxa{j}(1:5)];
disp(['Simulating spectral mixture for combination (' num2str(i) ',' ...
      num2str(j) '): ' mixCode{i,j}])
[SimSpec(i,j).outputRoot,SimSpec(i,j).fractions,...
 SimSpec(i,j).inputMix,SimSpec(i,j).fracDiff,SimSpec(i,j).rmse,...
 SimSpec(i,j).meanFracError,SimSpec(i,j).stdFracError,...
 SimSpec(i,j).meanRmse,SimSpec(i,j).stdRmse] ...
            =   HABsimSpec(libWater,wvl,Taxa,iTaxa2mix1,iTaxa2mix2,...
                     iWater,mixCode{i,j},dataDir,minEmFrac,maxEmFrac,...
                        minShadeEmFrac,maxShadeEmFrac,maxRMSE);
set(gcf,'name','AnaMicroSimMixExample')
saveas(gcf,fullfile('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SimSpec',get(gcf,'name')),'fig')

%%
%% Loop over all taxa combinations, create simulated mixtures, perform MESMA, and summarize output
for i = 1:numEm
    for j = 1:numEm
        if i < j
            iTaxa2mix1  =   i;
            iTaxa2mix2  =   j;
            mixCode{i,j}=   [Taxa{i}(1:5) Taxa{j}(1:5)];
            disp(['Simulating spectral mixture for combination (' num2str(i) ',' ...
                  num2str(j) '): ' mixCode{i,j}])
            [SimSpec(i,j).outputRoot,SimSpec(i,j).fractions,...
             SimSpec(i,j).inputMix,SimSpec(i,j).fracDiff,SimSpec(i,j).rmse,...
             SimSpec(i,j).meanFracError,SimSpec(i,j).stdFracError,...
             SimSpec(i,j).meanRmse,SimSpec(i,j).stdRmse] ...
                        =   HABsimSpec(libWater,wvl,Taxa,iTaxa2mix1,iTaxa2mix2,...
                                 iWater,mixCode{i,j},dataDir,minEmFrac,maxEmFrac,...
                                    minShadeEmFrac,maxShadeEmFrac,maxRMSE);
            meanFracErr(i,j,:)  =   SimSpec(i,j).meanFracError;
            stdFracErr(i,j,:)   =   SimSpec(i,j).stdFracError;
            meanRmse(i,j)       =   SimSpec(i,j).meanRmse;
            stdRmse(i,j)        =   SimSpec(i,j).stdRmse;
            close all
        end
    end
end

%%
%% Save what we have so far
clear i j ans iTaxa* 
save HABsimSpecAllPairs.mat
% Delete all intermediate simulated mixture images and MESMA files in Windows
% explorer to save disk space

%%
%% Make abbreviated fraction names
% TaxaShort  =   Taxa;
% for i = 1:length(TaxaShort)
%     TaxaShort{i} =    TaxaShort{i}(1:4);
% end
TaxaShort  =   TaxaLegend;
for i = 1:length(TaxaShort)-1
    TaxaShort{i} =    [TaxaShort{i}(1:9) '}}'];
end


%%
%% Visualize results: Mean fraction errors
% For now, average over the three end members fraction errors to get one value
% for each mixture
% Make a mask to only show the actual combinations
mask    =   zeros(numEm,numEm);
for i = 1:numEm
    for j = 1:numEm
        if i < j
            mask(i,j) = 1;
        end
    end
end
mask    =   logical(mask);
figure
t = tiledlayout(2,2);
nexttile
imagesc(mean(meanFracErr,3),'alphadata',mask); 
axis equal; axis tight
crameri('imola'); colorbar
set(gca,'Xtick',1:numEm,'XtickLabel',TaxaShort(1:numEm))
set(gca,'Ytick',1:numEm,'YtickLabel',TaxaShort(1:numEm))
title('(a) Mean fraction errors')
% set(gcf,'name','SimSpecMeanFracErr')
% saveas(gcf,get(gcf,'name'),'fig')
% saveas(gcf,get(gcf,'name'),'jpg')

%%
%% Visualize results: standard deviation of fraction errors
% For now, average over the three end members fraction errors to get one value
% for each mixture
% figure
nexttile
imagesc(mean(stdFracErr,3),'alphadata',mask); 
axis equal; axis tight
crameri('imola'); colorbar
set(gca,'Xtick',1:numEm,'XtickLabel',TaxaShort(1:numEm))
set(gca,'Ytick',1:numEm,'YtickLabel',TaxaShort(1:numEm))
title('(b) Std. dev. of fraction errors')
% set(gcf,'name','SimSpecStdFracErr')
% saveas(gcf,get(gcf,'name'),'fig')
% saveas(gcf,get(gcf,'name'),'jpg')

%%
%% Visualize results: mean MESMA RMSE
% figure
nexttile
imagesc(meanRmse,'alphadata',mask); 
axis equal; axis tight
crameri('imola'); colorbar
set(gca,'Xtick',1:numEm,'XtickLabel',TaxaShort(1:numEm))
set(gca,'Ytick',1:numEm,'YtickLabel',TaxaShort(1:numEm))
title('(c) Mean RMSE')
% set(gcf,'name','SimSpecMeanRMSE')
% saveas(gcf,get(gcf,'name'),'fig')
% saveas(gcf,get(gcf,'name'),'jpg')


%%
%% Visualize results: standard deviation of MESMA RMSE
% figure
nexttile
imagesc(stdRmse,'alphadata',mask); 
axis equal; axis tight
crameri('imola'); colorbar
set(gca,'Xtick',1:numEm,'XtickLabel',TaxaShort(1:numEm))
set(gca,'Ytick',1:numEm,'YtickLabel',TaxaShort(1:numEm))
title('(d) Std. dev. of RMSE')
% set(gcf,'name','SimSpecStdRMSE')
% saveas(gcf,get(gcf,'name'),'fig')
% saveas(gcf,get(gcf,'name'),'jpg')

%%
%% Save the figure to disk
set(gcf,'name','SimSpecErrorMatrices')
saveas(gcf,get(gcf,'name'),'fig')
saveas(gcf,get(gcf,'name'),'jpg')










%%
%% **** 10/5/2021: ADJUST RMSE TO AVOID FALSE POSITIVES WITHOUT LOSING SENSITIVITY
% We want to avoid classifying taxa that are not actually present just because
% the MESMA algorithm wants to classify something.  One way of preventing this
% kind of false negative might be to make the RMSE threshold more stringent,
% smaller than the default 0.025 we've used so far.  However, if we make the
% RMSE cutoff too small we might also risk losing the sensitivity to detect taxa
% that are actually present and making too many unclassified pixels. So, there
% is a balance to be struck and we need to essentially calibrate an appropriate
% RMSE threshold.  To do so, let's select a water body where we know that the
% algal taxa sampled in the field was detected in the image (AFA in UKL) and
% another water body where the taxa sampled in the field is not included in our
% library and so we should get unclassified pixels (Detroit Lake). We can apply
% the same RMSE to both images and hopefully retain the correct classifications
% in UKL while avoiding misleading classifications in Detroit Lake.
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')

%%
%% Load previous results from DL20200813 and repeat with a more stringent RMSE
% load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813\DL20200813.mat')
% Initially we used a max RMSE threshold of 0.025, so let's reduce this to
% something very small like 0.01 and re-run the MESMA
maxRMSE     =   0.002;
% Keep everything else the same and re-run the MESMA
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
%%                      
%% DL20200813 low RMSE post-SMASH
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   


%%
%% Now check to see if this low RMSE threshold causes us to lose sensitivity at UKL
% close all
% clear all
% % Load previous results from UK20200810 and repeat with a more stringent RMSE
% load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200828.mat')
% Initially we used a max RMSE threshold of 0.025, so let's reduce this to
% something very small like 0.01 and re-run the MESMA
maxRMSE     =   0.025;
% Noticed that date code was wrong, so fix here:
dateCode    =   '20200810';
% Keep everything else the same and re-run the MESMA
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
%%                      
%% UK20200810 low RMSE post-SMASH
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);

%%
%% Looks like we're not going to have a universal RMSE threshold, so ...
% An RMSE cutoff that was stringent enough to avoid spurious classifications at
% Detroit was too stringent to allow for any classification, let alone a correct
% classification, at Upper Klamath. So ...

% We will need to write some new code to repeat the MESMA for different RMSE
% cutoffs for a different lake and see where we start to get correct
% classifications at Klamath again and where we start to get spurious
% classifications at Detroit

%%
%% Run new function tuneRmse.m for Detroit Lake test case
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813\DL20200813.mat')
% % Initially we used a max RMSE threshold of 0.025, so let's reduce this to
% % something very small like 0.01 and re-run the MESMA
maxRmse     =   0.001:0.0005:0.01;
runDelay    =   30;
[Taxa,MesmaOut,classProp] = ...
    tuneRmse(maxRmse,runDelay,dataDir,siteCode,dateCode,libraryFile,...
             classNameField,imgSubFile,shadeFile,imgScaleFactor,libScaleFactor,...
             minEmFrac,maxEmFrac,minShadeEmFrac,maxShadeEmFrac,waterMask,R);

%%
%% Save what we have so far
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813\DL20200813rmseSensitivity.mat','-v7.3')

%%
%% Now for UKL, too, which might have the opposite problem
% load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200828.mat')
maxRmse     =   0.001:0.001:0.025;
runDelay    =   180;
[Taxa,MesmaOut,classProp] = ...
    tuneRmse(maxRmse,runDelay,dataDir,siteCode,dateCode,libraryFile,...
             classNameField,imgSubFile,shadeFile,imgScaleFactor,libScaleFactor,...
             minEmFrac,maxEmFrac,minShadeEmFrac,maxShadeEmFrac,waterMask,R);

%%
%% Save what we have so far
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200828rmseSensitivity.mat','-v7.3')

%%
%% NOTES FOR NEXT TIME
% Need to extend this analysis to a broader range of RMSE thresholds and for
% multiple sites to see if we can hone in on some universal value that works
% well in all cases


%% **** 10/11/2021: FIRST CLEAN UP MESMA OUTPUT FILES FROM A DIRECTORY
outDir      =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810';
oldFiles    =   dir([outDir '\*Mesma*']);
for i = 1:length(oldFiles)
    delete(fullfile(outDir,oldFiles(i).name))
end

%%
%% First let's try the other two dates at Detroit, starting with DL20190819
% load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819\DL20190819.mat')
% Initially we used a max RMSE threshold of 0.025, so let's reduce this to
% something very small like 0.01 and re-run the MESMA
maxRmse     =   0.001:0.0005:0.01;
runDelay    =   20;
[Taxa,MesmaOut,classProp] = ...
    tuneRmse(maxRmse,runDelay,dataDir,siteCode,dateCode,libraryFile,...
             classNameField,imgSubFile,shadeFile,imgScaleFactor,libScaleFactor,...
             minEmFrac,maxEmFrac,minShadeEmFrac,maxShadeEmFrac,waterMask,R);

%%
%% Save what we have so far and clean up MESMA output files
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819\DL20190819rmseSensitivity.mat','-v7.3')
outDir      =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819';
oldFiles    =   dir([outDir '\*Mesma*']);
for i = 1:length(oldFiles)
    delete(fullfile(outDir,oldFiles(i).name))
end

%%
%% Now for DL20200828
% load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828\DL20200828.mat')
% Initially we used a max RMSE threshold of 0.025, so let's reduce this to
% something very small like 0.01 and re-run the MESMA
maxRmse     =   0.001:0.0005:0.01;
runDelay    =   20;
[Taxa,MesmaOut,classProp] = ...
    tuneRmse(maxRmse,runDelay,dataDir,siteCode,dateCode,libraryFile,...
             classNameField,imgSubFile,shadeFile,imgScaleFactor,libScaleFactor,...
             minEmFrac,maxEmFrac,minShadeEmFrac,maxShadeEmFrac,waterMask,R);

%%
%% Save what we have so far and clean up MESMA output files
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828\DL20200828rmseSensitivity.mat','-v7.3')
outDir      =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828';
oldFiles    =   dir([outDir '\*Mesma*']);
for i = 1:length(oldFiles)
    delete(fullfile(outDir,oldFiles(i).name))
end

%%
%% Back East to OWASCO - OW20190805
% load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805\OW20190805.mat')
% Initially we used a max RMSE threshold of 0.025, so let's reduce this to
% something very small like 0.01 and re-run the MESMA
maxRmse     =   0.001:0.0005:0.01;
runDelay    =   30;
[Taxa,MesmaOut,classProp] = ...
    tuneRmse(maxRmse,runDelay,dataDir,siteCode,dateCode,libraryFile,...
             classNameField,imgSubFile,shadeFile,imgScaleFactor,libScaleFactor,...
             minEmFrac,maxEmFrac,minShadeEmFrac,maxShadeEmFrac,waterMask,R);

%%
%% Save what we have so far and clean up MESMA output files
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805\OW20190805rmseSensitivity.mat','-v7.3')
outDir      =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805';
oldFiles    =   dir([outDir '\*Mesma*']);
for i = 1:length(oldFiles)
    delete(fullfile(outDir,oldFiles(i).name))
end

%%
%% Second date for OWASCO - OW20200823
% load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823\OW20200823.mat')
% Initially we used a max RMSE threshold of 0.025, so let's reduce this to
% something very small like 0.01 and re-run the MESMA
maxRmse     =   0.001:0.0005:0.01;
runDelay    =   45;
[Taxa,MesmaOut,classProp] = ...
    tuneRmse(maxRmse,runDelay,dataDir,siteCode,dateCode,libraryFile,...
             classNameField,imgSubFile,shadeFile,imgScaleFactor,libScaleFactor,...
             minEmFrac,maxEmFrac,minShadeEmFrac,maxShadeEmFrac,waterMask,R);

%%
%% Save what we have so far and clean up MESMA output files
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823\OW20200823rmseSensitivity.mat','-v7.3')
outDir      =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823';
oldFiles    =   dir([outDir '\*Mesma*']);
for i = 1:length(oldFiles)
    delete(fullfile(outDir,oldFiles(i).name))
end

%%
%% **** 10/11/2021: MOSAIC GRAPEVINE
% Done in ENVI to create a new image file 
    % C:\Users\cjl\OneDrive -
    % DOI\HABs\SMASH\Grapevine\GV20200813mosaic\rawImage\GV20200813mosaic.tif
% Need to do an initial SMASH run to set up the inputs we'll need for tuneRmse.m    

%%
%% GV20200813mosaic pre-SMASH
% Load generic inputs from initial run before we mosaicked 
% load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813\GV20200813.mat')
close all
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813mosaic';
siteCode    =   'GV';
dateCode    =   '20200813';
siteDateCode=   [siteCode dateCode];
% WGS 84 UTM Zone 14S is EPSG 32614;
epsgCode    =   32614;
imgFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813mosaic\rawImage\GV20200813mosaic.tif';
wvl2keep    =   [400 805];
maskInput   =   [];
% Used [0 2100] for threshold-based mask developed from NDWI image
[R,ndwi,waterMask,ndci,ci,wvlSub,imgSpecSub,imgSubFile] = ...
      preSMASH(dataDir,siteCode,dateCode,epsgCode,wvl2keep,imgFile,maskInput); %#ok<*ASGLU>
imtool close all
close(1)
figDir      =   [dataDir '\figs'];
if ~isfolder(figDir); mkdir(figDir); end
figure(2)
set(gcf,'name',[siteDateCode '_NDCI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_CI'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

%%
%% GV20200813mosaic run SMASH - most of these inputs won't change from site to site
% This will change from scene to scene but can be updated programmatically
% Name of image file
imgSubFile      =   fullfile(dataDir,[siteDateCode 'specSub.tif']);

% The rest should be consistent for all sites
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\LibrarySmoothed.sli';
% Metadata field with end-member names
classNameField  =   'Class';
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
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\WaterSmoothed.sli';

% Call the function
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);
                      
%%                      
%% GV20200813mosaic post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))

%% 
%% Now run tuneRmse for Grapevine
% Initially we used a max RMSE threshold of 0.025, so let's reduce this to
% something very small like 0.01 and re-run the MESMA
maxRmse     =   0.001:0.0005:0.01;
runDelay    =   30;
[Taxa,MesmaOut,classProp] = ...
    tuneRmse(maxRmse,runDelay,dataDir,siteCode,dateCode,libraryFile,...
             classNameField,imgSubFile,shadeFile,imgScaleFactor,libScaleFactor,...
             minEmFrac,maxEmFrac,minShadeEmFrac,maxShadeEmFrac,waterMask,R);

%%
%% Save what we have so far and clean up MESMA output files
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813mosaic\GV20200813rmseSensitivity.mat','-v7.3')
outDir      =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813mosaic';
oldFiles    =   dir([outDir '\*Mesma*']);
for i = 1:length(oldFiles)
    delete(fullfile(outDir,oldFiles(i).name))
end

%%
%% NOTES FOR NEXT TIME
% Think about tuneRmse results for all sites/dates and try to identify a commmon
% max RMSE threshold we can apply to all data sets
% Use this max RMSE to repeat a set of runs with common inputs of all parameters
% Extract fractions from field sample locations
% Maybe work on a more general legend placement

%%
%% **** 10/12/2021: SUMMARIZE FIELD-SAMPLED TAXA AND APPROPRIATE MAX RMSE CUTOFF
% OW20190805: Synecocchus not in library, so unclassfied is correct --> 0.0025

% OW20200823: Microsystis --> 0.01

% GV20200813: Raphidiopsis (aka Cylindrospermopsis) not in library, unless this
% is the same as Cylindrospermum; both are in the family Aphanizomenonaceae, so
% let's use 0.01 which will give us mostly AFA but also about 20% unclassified

% DL20190819: Diatoms not in library, so unclassified is correct --> 0.002

% DL20200813: Dolichospermum not in library, so unclassified is corrrect --> 0.002

% DL20200828: Gloeotrichia in blowout creek arm only, but image infers
% Oscillatoria, so mostly unclassified is probably right --> 0.004

% UK20200828: AFA in field, but requires a higher max RMSE to detect in image
% --> 0.02

%%
%% OW20190805: repeat run with max RMSE of 0.0025
% Load previous results so we have the outputs from preSMASH
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805\OW20190805.mat')
% Now repeat with the new maxRMSE
maxRMSE     =   0.0025;
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);

%% OW20190805: post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(9)
set(gcf,'name',[siteDateCode '_RMSEbyTaxa'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');


% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))

%%
%% OW20200823: repeat run with max RMSE of 0.01
% Load previous results so we have the outputs from preSMASH
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823\OW20200823.mat')
% Now repeat with the new maxRMSE
maxRMSE     =   0.01;
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);

%% OW20200823: post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(9)
set(gcf,'name',[siteDateCode '_RMSEbyTaxa'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))

%%
%% GV20200813: repeat run with max RMSE of 0.01
% Load previous results so we have the outputs from preSMASH
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813mosaic\GV20200813.mat')
% Now repeat with the new maxRMSE
maxRMSE     =   0.01;
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);

%% GV20200813: post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(9)
set(gcf,'name',[siteDateCode '_RMSEbyTaxa'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))

%%

%% DL20190819: repeat run with max RMSE of 0.002
% Load previous results so we have the outputs from preSMASH
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819\DL20190819.mat')
% Now repeat with the new maxRMSE
maxRMSE     =   0.002;
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);

%% DL20190819: post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(9)
set(gcf,'name',[siteDateCode '_RMSEbyTaxa'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))

%%

%% DL20200813: repeat run with max RMSE of 0.002
% Load previous results so we have the outputs from preSMASH
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813\DL20200813.mat')
% Now repeat with the new maxRMSE
maxRMSE     =   0.002;
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);

%% DL20200813: post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(9)
set(gcf,'name',[siteDateCode '_RMSEbyTaxa'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))

%%

%% DL20200828: repeat run with max RMSE of 0.004
% Load previous results so we have the outputs from preSMASH
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828\DL20200828.mat')
% Now repeat with the new maxRMSE
maxRMSE     =   0.004;
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);

%% DL20200828: post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(9)
set(gcf,'name',[siteDateCode '_RMSEbyTaxa'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))

%%
%% UK20200810: repeat run with max RMSE of 0.02
% Load previous results so we have the outputs from preSMASH
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810.mat')
% Now repeat with the new maxRMSE
maxRMSE     =   0.02;
dateCode   =    '20200810';
imgSubFile =    'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810specSub.tif';
outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                          classNameField,imgSubFile,shadeFile,...
                          imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                          minShadeEmFrac,maxShadeEmFrac,maxRMSE);

%%
%% UK20200810: post-SMASH - save figures and a .mat file with results
close all
[Taxa,model,modelSum,taxaHist,taxaId,classMask,fractions,taxaIdR,...
     taxaIdG,rgbFrac,rmse,meanRmse] = postSMASH(imgSubFile,wvlSub,R,siteCode,...
                                       dateCode,waterMask,libraryFile,outputRoot);                                   
figure(1)
set(gcf,'name',[siteDateCode '_RGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(2)
set(gcf,'name',[siteDateCode '_MESMAclassification'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(3)
set(gcf,'name',[siteDateCode '_MESMAclassHistogram'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(4)
set(gcf,'name',[siteDateCode '_selectedTaxaMap'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(5)
set(gcf,'name',[siteDateCode '_selectedTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(6)
set(gcf,'name',[siteDateCode '_allTaxaFractions'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(7)
set(gcf,'name',[siteDateCode '_taxaFractionRGB'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(8)
set(gcf,'name',[siteDateCode '_RMSE'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');
figure(9)
set(gcf,'name',[siteDateCode '_RMSEbyTaxa'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig');
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg');

% Also save all results to a MAT file
save(fullfile(dataDir,[siteDateCode '.mat']))

%%






%% **** 10/13/2021: EXTRACT MESMA FRACTIONS FROM FIELD SAMPLE LOCATIONS
% Link to spreadsheet with matchups:
    % https://doimspp.sharepoint.com/:x:/r/sites/GS-WMA-RSWQ/_layouts/15/Doc.aspx?sourcedoc=%7BD1CB64AA-5153-4FCD-93EC-F3F6EC8B1054%7D&file=HS_satellite_in-situ_match_ups.xlsx&action=default&mobileredirect=true&cid=f68501ea-0db6-413d-a885-a595caa63983
% Also saved locally as:
    % C:\Users\cjl\OneDrive - DOI\HABs\SMASH\HS_satellite_in-situ_match_ups.xlsx
% Use the function we developed for this purpose: extractFractions.m

%%
%% OW20190805
% Field-sampled taxa and percent biovolume: Synechococcus (coccoid cyano, <1 um), 63%
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805\OW20190805.mat');
fieldLat    =   42.89090194;
fieldLon    =   -76.52675986;
disp(['ESPG code = ' num2str(epsgCode)])
fracAtField =   extractFractions(fieldLat,fieldLon,siteCode,dateCode,...
                                 epsgCode,fractions,Taxa,R);
set(gcf,'name',[siteCode dateCode 'fracAtField'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg')
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805\OW20190805.mat');

%%
%% OW20200823
% Field-sampled taxa and percent biovolume: Microcystis spp (M. aeriginosa, M. flos-aquae), 88 %
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823\OW20200823.mat');
fieldLat    =   42.89090194;
fieldLon    =   -76.52675986;
disp(['ESPG code = ' num2str(epsgCode)])
fracAtField =   extractFractions(fieldLat,fieldLon,siteCode,dateCode,...
                                 epsgCode,fractions,Taxa,R);
set(gcf,'name',[siteCode dateCode 'fracAtField'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg')
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823\OW20200823.mat');

%%
%% GV20200813
% Field-sampled taxa and percent biovolume: Raphidiopsis (aka Cylindrospermopsis), 24%
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813mosaic\GV20200813.mat');
fieldLat    =   33.0132099;
fieldLon    =   -97.15011136;
disp(['ESPG code = ' num2str(epsgCode)])
fracAtField =   extractFractions(fieldLat,fieldLon,siteCode,dateCode,...
                                 epsgCode,fractions,Taxa,R);
set(gcf,'name',[siteCode dateCode 'fracAtField'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg')
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813mosaic\GV20200813.mat');

%%
%% DL20190819
% Field-sampled taxa and percent biovolume: Diatoms (Asterionella, Fragilaria crotonensis), 50-55 %
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819\DL20190819.mat');
fieldLat    =   44.71913718;
fieldLon    =   -122.2460034;
disp(['ESPG code = ' num2str(epsgCode)])
fracAtField =   extractFractions(fieldLat,fieldLon,siteCode,dateCode,...
                                 epsgCode,fractions,Taxa,R);
set(gcf,'name',[siteCode dateCode 'fracAtField'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg')
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819\DL20190819.mat');

%%
%% DL20200813
% Field-sampled taxa and percent biovolume: Dolichospermum, 52 %
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813\DL20200813.mat');
fieldLat    =   44.720215047318206;
fieldLon    =    -122.24836558729321;
disp(['ESPG code = ' num2str(epsgCode)])
fracAtField =   extractFractions(fieldLat,fieldLon,siteCode,dateCode,...
                                 epsgCode,fractions,Taxa,R);
set(gcf,'name',[siteCode dateCode 'fracAtField'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg')
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813\DL20200813.mat');

%%
%% DL20200828
% Field-sampled taxa and percent biovolume: Gloeotrichia, 83 %
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828\DL20200828.mat');
fieldLat    =   44.67546;
fieldLon    =   -122.185;
disp(['ESPG code = ' num2str(epsgCode)])
fracAtField =   extractFractions(fieldLat,fieldLon,siteCode,dateCode,...
                                 epsgCode,fractions,Taxa,R);
set(gcf,'name',[siteCode dateCode 'fracAtField'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg')
save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828\DL20200828.mat');

%%
%% UK20200810 - two sites in this lake
% Field-sampled taxa and percent biovolume: AFA field verified/lab hyperspectral AFA
load('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810.mat');
dateCode    =   '20200810'; 
Taxa        =   ['Unclassified'; Taxa];
fieldLatRSP =   42.344214;
fieldLonRSP =   -121.858296;
fieldLatWRM =   42.463547;
fieldLonWRM =   -121.95992;
disp(['ESPG code = ' num2str(epsgCode)])
fracAtFieldRSP  =   extractFractions(fieldLatRSP,fieldLonRSP,siteCode,dateCode,...
                                 epsgCode,fractions,Taxa,R);
set(gcf,'name',[siteCode dateCode 'fracAtFieldRSP'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg')

fracAtFieldWRM  =   extractFractions(fieldLatWRM,fieldLonWRM,siteCode,dateCode,...
                                 epsgCode,fractions,Taxa,R);
set(gcf,'name',[siteCode dateCode 'fracAtFieldWRM'])
saveas(gcf,fullfile(figDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(figDir,get(gcf,'name')),'jpg')

save('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810.mat');

%%



%% **** 10/18/2021: UPDATE FIGURES FOR EACH SITE/DATE FOR MANUSCRIPT
% Write a new function to read the existing figures, make a few changes, and
% save them as new .fig, .jpg, and .eps files. We can compile the .eps files in
% a series of LaTeX documents, one for each site/date, to include as supporting
% information for the manuscript. This approach will be cleaner and less labor
% intensive than the traditional PowerPoint figure dump.  Place the new figures
% into a separate "figsUpdate" folder within each site/date folder we already
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
%% **** 10/20/2021: CREATE FIGURES FOR MANUSCRIPT USING NEW updateFigs.m FUNCTION
%%
%% UK20200810: done as a prototype in developing updateFigs.m
matFile         =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810.mat';
scaleBarDist    =   10000;
scaleOffset     =   [5000 -5000];
cBarOffset      =   0.05;
legLabelOffset  =   0.02;
legBoxWidth     =   0.05;
updateFigs(matFile,scaleBarDist,scaleOffset,cBarOffset,legLabelOffset,legBoxWidth);

%%
%% DL20190819
matFile         =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20190819\DL20190819.mat';
scaleBarDist    =   5000;
scaleOffset     =   [0 -500];
cBarOffset      =   -0.03;
legLabelOffset  =   0.02;
legBoxWidth     =   0.05;
updateFigs(matFile,scaleBarDist,scaleOffset,cBarOffset,legLabelOffset,legBoxWidth);

%%
%% DL20200813
matFile         =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200813\DL20200813.mat';
scaleBarDist    =   5000;
scaleOffset     =   [-100 -600];
cBarOffset      =   0.03;
legLabelOffset  =   0.02;
legBoxWidth     =   0.05;
updateFigs(matFile,scaleBarDist,scaleOffset,cBarOffset,legLabelOffset,legBoxWidth);

%%
%% DL20200828
matFile         =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Detroit\DL20200828\DL20200828.mat';
scaleBarDist    =   5000;
scaleOffset     =   [-100 -600];
cBarOffset      =   0.03;
legLabelOffset  =   0.02;
legBoxWidth     =   0.05;
updateFigs(matFile,scaleBarDist,scaleOffset,cBarOffset,legLabelOffset,legBoxWidth);

%%
%% GV20200813
matFile         =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Grapevine\GV20200813mosaic\GV20200813.mat';
scaleBarDist    =   5000;
scaleOffset     =   [-100 -600];
cBarOffset      =   0.07;
legLabelOffset  =   0.02;
legBoxWidth     =   0.05;
updateFigs(matFile,scaleBarDist,scaleOffset,cBarOffset,legLabelOffset,legBoxWidth);

%%
%% OW20190805
matFile         =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20190805\OW20190805.mat';
scaleBarDist    =   1000;
scaleOffset     =   [0 -200];
cBarOffset      =   0.07;
legLabelOffset  =   0.02;
legBoxWidth     =   0.05;
updateFigs(matFile,scaleBarDist,scaleOffset,cBarOffset,legLabelOffset,legBoxWidth);

%%
%% OW20200823
matFile         =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823\OW20200823.mat';
scaleBarDist    =   5000;
scaleOffset     =   [0 -200];
cBarOffset      =   0.07;
legLabelOffset  =   0.02;
legBoxWidth     =   0.05;
updateFigs(matFile,scaleBarDist,scaleOffset,cBarOffset,legLabelOffset,legBoxWidth);

%%
%% **** 10/20/2021: CREATE UKL FIGURES TO INCLUDE IN MANUSCRIPT ITSELF
load("C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810.mat")
newFigDir   =   "C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\updateFigs";

%%
%% Three panel figure with RGB, NDCI, and CI
set(0,'defaultFigureWindowStyle','normal')
hFig        =   figure('Position',[50 50 1600 800]);
% Set scale bar length and shift relative to default position
scaleBarDist=   10000;
scaleOffset =   [5000 -5000];

% RGB image
hAx1        =       subplot(1,3,1);
% imgSubFile  =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\UK20200810specSub.tif';
% Hcube       =   hypercube(imgSubFile,wvlSub);
% [rgb,iRgb]  =   colorize(Hcube,'Method','rgb','ContrastStretching',1);
% rgbWvl      =   wvlSub(iRgb);
% % Flip to account for georeferencing
% rgb         =   flipud(rgb);
imagesc(rgb,'XData',R.XWorldLimits,'YData',R.YWorldLimits,'alphadata',waterMask);
axis image; axis xy
% title({[siteDateCode ' contrast-stretched RGB image'],...
%        ['R: ' num2str(rgbWvl(1)) ' nm, G: ' num2str(rgbWvl(2)) ' nm, B: ' ...
%        num2str(rgbWvl(3)) ' nm']});
title('(a) RGB image')
[hLine,hText]       =   scaleBar(scaleBarDist);
hLine.XData         =   hLine.XData + scaleOffset(1);
hLine.YData         =   hLine.YData + scaleOffset(2);
hText.Position(1)   =   hText.Position(1) + scaleOffset(1);
hText.Position(2)   =   hText.Position(2) + scaleOffset(2);
hText.String        =   [num2str(scaleBarDist/1000) ' km'];
axis off

% NDCI image
hAx2                =    subplot(1,3,2);
imagesc(ndci,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
    'AlphaData',waterMask);
axis image; axis xy;
crameri('-bamako')
hCbar1              =    colorbar;
hCbar1.Position     =    [0.62    0.3044    0.0133    0.4086];
% title([siteDateCode ': Normalized Difference Chlorophyll Index'])
title('(b) Normalized Difference Chlorophyll Index (NDCI)')
axis off

% CI image
hAx3                =   subplot(1,3,3);
imagesc(ci,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
    'AlphaData',waterMask);
axis image; axis xy; axis off
crameri('-bamako')
hCbar2              =    colorbar;
hCbar2.Position     =    [0.9    0.3044    0.0133    0.4086];
% title([siteDateCode ': Cyanobacterial Index'])
title('(c) Cyanobacterial Index (CI)')

set(gcf,'name',[siteDateCode '_RGB_NDCI_CI'])
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'fig')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'jpg')
saveas(gcf,fullfile(newFigDir,get(gcf,'name')),'epsc2')

%% Output from tuneRmse.m - can use existing figure
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\updateFigs\UK20200810_maxRMSEsensitivity.fig and .eps

%% Use existing MESMA-based classification and histogram
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\updateFigs\UK20200810_MESMAclassification.fig and .eps
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\updateFigs\UK20200810_MESMAclassHistogram.fig and .eps

%% Use existing fraction image mosaic for all endmembers
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\updateFigs\UK20200810_allTaxaFractions.fig and .eps

%% Use existing endmember fraction composite
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\updateFigs\UK20200810_taxaFractionRGB.fig and .eps
% Include former title in caption

%% Use existing RMSE image and bar chart by genus
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\updateFigs\UK20200810_RMSE.fig and .eps
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\UKL\UK20200810\updateFigs\UK20200810_RMSEbyTaxa.fig and .eps


%% Testing out difftool in git