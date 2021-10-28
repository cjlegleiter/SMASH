%% Assess MESMA feasibility by quantifying spectral contrast and simulating mixtures
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 06/09/2021
% See some prototyping code in: C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SMASHowasco.m
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SimSpec\HABsimSpecLog.m


%% Begin by computing metrics of spectral contrast for end members in algal library
% The MATLAB Hyperspectral Imaging Toolbox has several spectra for this purpose
%   web(fullfile(docroot, 'images/hyperspectral-image-processing.html?s_tid=CRUX_lftnav'))

% Import the library
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
hLabel  =   uilistbox(hFig,'Items',TaxaNums,'Position',[5 5 290 290]);
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
title('Merged algae and water spectral library')
legend(TaxaNums)


%% Metrics of separability in a numEm X numEm matrix
% The Hyperspectral Imaging Toolbox has several measures of spectral similarity
% available but the best summary might be the normalized spectral similiarity
% score (NS3), which combines both the overall brightness (Euclidean distance)
% and spectral shape (spectral angle) in a summary score. See documentation for
% the function ns3: "The test spectrum with the minimum score value indicates
% highest similarity to the reference endmember. On the other hand, the test
% spectrum with the maximum score value has the highest spectral variability."
numEm   =   size(library,1);
ns3score=   zeros(numEm,numEm);
for i = 1:numEm
    for j = 1:numEm
        if i < j
            ns3score(i,j)   =   ns3(library(i,:),library(j,:));
        end
    end
end
figure
imagesc(ns3score,'alphadata',logical(ns3score)); 
axis equal; axis tight
crameri('imola'); colorbar
set(gca,'Xtick',1:14,'XtickLabel',Taxa)
set(gca,'Ytick',1:14,'YtickLabel',Taxa)
title('Normalized spectral similarity scores for algal taxa and water')
% Looks like anabena is relatively distinct from most of the other taxa by this
% metric, so let's try making a mixture of anabena vs. microcystis, which is
% actually harmful and has been observed in Lake Owasco
clear i j ans h*


%% **** Create a box of simulated spectral mixtures with two algal taxa and water


%% Specify taxa to mix and plot an example mixed spectrum
iTaxa2mix1  =   1;  % Anabena
iTaxa2mix2  =   9;  % Microcystis
iWater      =   14;
fracTaxa1   =   0.5;
fracTaxa2   =   0.4;
fracWater   =   1-fracTaxa1-fracTaxa2;
mix         =   library(iTaxa2mix1,:)*fracTaxa1 + ...
                library(iTaxa2mix2,:)*fracTaxa2 + ...
                library(iWater,:)*fracWater;
% Make the library figure active            
hold on
plot(wvlSub,mix,'k--','linewidth',2)
% Looks like this is working OK


%% Set up a hypotehtical fraction images for known simulated mixtures
% 100 rows each for a water fraction of 0 to 1 in steps of 0.01
fracWater   =   (0:0.01:1)';
% Fractions for the two algae evenly spaced between the two end members for the
% non-water portion of each pixel
fracTaxa1   =   0:0.01:1;
fracTaxa2   =   1-fracTaxa1;
% Now allocate the water fractions to an array
fracWaterImg=   zeros(length(fracWater),length(fracWater));
for i = 1:length(fracWater)
    fracWaterImg(i,:)   =   fracWater(i);
end
figure
subplot(1,3,1)
imshow(fracWaterImg)
colorbar
title('Water Fraction')
% Same kind of thing for the taxa fractions, adjusted to account for the water
% fraction already allocated for that row
fracTaxa1img=   zeros(length(fracWater),length(fracWater));
fracTaxa2img=   zeros(length(fracWater),length(fracWater));
for i = 1:length(fracWater)
    for j = 1:length(fracTaxa1)
        fracTaxa1img(i,j)   =   fracTaxa1(j)*(1-fracWater(i));
        fracTaxa2img(i,j)   =   fracTaxa2(j)*(1-fracWater(i));    
    end
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
simSpec     =   zeros(length(fracWater),length(fracWater),length(wvlSub));
for i = 1:length(fracWater)
   for j = 1:length(fracWater)
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
figure
simCube     =   hypercube(simSpec,wvlSub);
simRgb      =   colorize(simCube,"Method","rgb","ContrastStretching",false);
imagesc(simRgb,'Xdata',fracWater,'Ydata',fracTaxa1);
axis image;
xlabel({'Taxa 1 fraction','Taxa 2 fraction = 1 - Taxa 1 fraction'})
ylabel('Water fraction')
title({'RGB image for simulated spectral mixtures:',...
       ['Taxa 1 = ' Taxa{iTaxa2mix1} ', Taxa 2 = ' Taxa{iTaxa2mix2}]})


%% Next write this simulated image out as a TIFF we can use for MESMA
siteDateCode=   'AnaMicro';
dataDir     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SimSpec';
simFile     =   fullfile(dataDir,[siteDateCode 'SimSpec.tif']);
% Make up some arbitrary spatial referencing info for export as a GeoTIFF
R           =   maprefcells([-0.5 100.5],[-0.5 100.5],size(fracWaterImg),...
                            'ColumnsStartFrom','north'); % To avoid having to flip
% For Owasco, WGS 84 UTM Zone 18N is EPSG 32618
epsgCode    =   32618;
geotiffwrite(simFile,simSpec,R,'CoordRefSysCode',epsgCode);
% % MATLAB to envi as an alternative?
% simFile     =   fullfile(dataDir,[siteDateCode 'SimSpec.dat']);
% OutHeader   =   matlab2envi(uint16(simSpec*10000),simFile,[],...
%                             -0.5:100.5,-0.5:100.5,wvlSub);


%% Same thing for the fraction images
% Also export the fraction images for the two taxa and water
fracTax1file=   fullfile(dataDir,[siteDateCode 'FracTaxa1.tif']);
fracTax2file=   fullfile(dataDir,[siteDateCode 'FracTaxa2.tif']);
fracWatrFile=   fullfile(dataDir,[siteDateCode 'FracWater.tif']);
geotiffwrite(simFile,uint8(fracTaxa1img*100),R,'CoordRefSysCode',epsgCode);
geotiffwrite(simFile,uint8(fracTaxa2img*100),R,'CoordRefSysCode',epsgCode);
geotiffwrite(simFile,uint8(fracWaterImg*100),R,'CoordRefSysCode',epsgCode);


%% Subset the library to include only the three components in these mixtures
librarySub  =   library([iTaxa2mix1 iTaxa2mix2 iWater],:);
% Best way to proceed will be to just write a new ENVI header using the code in
% the new function specLibHeader.m after calling matlab2envi to actually write
% out the data
matlab2envi(librarySub,'AnaMicroWater.sli',outHeader,[],[],[]);
Header      =   specLibHeader(librarySub,'AnaMicroWater.sli',Taxa([1 9 14]),wvlSub);
% Also make the corresponding csv metadata file
Metadata    =   table(Taxa([iTaxa2mix1 iTaxa2mix2 iWater]),...
                      Taxa([iTaxa2mix1 iTaxa2mix2 iWater]),...
                      'VariableNames',{'Name','Class'});
writetable(Metadata,'AnaMicroWater.csv')


%% Subset the library to include only the two algal taxa in these mixtures
librarySub  =   library([iTaxa2mix1 iTaxa2mix2],:);
% Best way to proceed will be to just write a new ENVI header using the code in
% the new function specLibHeader.m after calling matlab2envi to actually write
% out the data
matlab2envi(librarySub,'AnaMicro.sli',[],[],[],wvlSub);
Header      =   specLibHeader(librarySub,'AnaMicro.sli',Taxa([iTaxa2mix1 iTaxa2mix2]),wvlSub); %#ok<*NASGU>
% Also make the corresponding csv metadata file
Metadata    =   table(Taxa([iTaxa2mix1 iTaxa2mix2]),...
                      Taxa([iTaxa2mix1 iTaxa2mix2]),...
                      'VariableNames',{'Name','Class'});
writetable(Metadata,'AnaMicro.csv')


%% Subset the library to include only the water in these mixtures
librarySub  =   library(iWater,:);
% Best way to proceed will be to just write a new ENVI header using the code in
% the new function specLibHeader.m after calling matlab2envi to actually write
% out the data
matlab2envi(librarySub,'Water.sli',[],[],[],wvlSub);
Header      =   specLibHeader(librarySub,'Water.sli',Taxa(iWater),wvlSub);
% Also make the corresponding csv metadata file
Metadata    =   table(Taxa(iWater),...
                      Taxa(iWater),...
                      'VariableNames',{'Name','Class'});
writetable(Metadata,'Water.csv')


%% Set parameters we will pass to the MESMA core algorithm
% Name of input library
libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SimSpec\AnaMicro.sli';
% Metadata field with end-member names
classNameField  =   'Class';
% Name of image file
imgSubFile      =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SimSpec\AnaMicroSimSpec.tif';
% Scale factors for library and image
libScaleFactor  =   1;
imgScaleFactor  =   1;
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
% % Name of shade end member file
shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SimSpec\Water.sli';


%% Programmatically generate new magic.bat file with actual call to mesma embedded
% Set up the actual call to the mesma routine itself first so we can embed it in
% the .bat file we're writing here
% Create a unique file name for the .bat file to be generated
timeStamp   =   datestr(now,'yyyymmdd_HHMMSS');
batFileName =   fullfile(dataDir,[siteDateCode '_CallMesma_' timeStamp '.bat']);
outputRoot  =   fullfile(dataDir,[siteDateCode '_MesmaOutput_' timeStamp]);
mesmaCmd    =   ['mesma "' libraryFile '" ' classNameField ' "' imgSubFile '"' ...
                 ' -l 3' ... Set complexity level to run all 3-EM models
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


%% Now call the new .bat file to run the mesma
% Specify the command as a string and then call it using the system function
command =   [batFileName ' &'];
status  =   system(command);
if status == 0
    disp('Issued programmatic call to MESMA, now awaiting results ...')
else
    disp('WARNING: MESMA did not launch as intended ...')
end


%% Looks promising if we use only the two algal taxa as EMs and then water as shade EM


%% Save what we have so far
clear ans fid i libSubCube outHeader OutHeader status t tmp
save HABsimSpec.mat


%% **** 6/17/2021: IMPORT MESMA RESULTS
% See prototyping code in SMASHowasco.m


%% *MESMA_fractions file contains the model's fractions, including shade
% The model’s fractions [nb of bands = nb of classes + 1], including a shade fraction
% Value of unmodeled pixels in output: 0
% Value of pixels with no data in output: 0
fractionsFile   =   [outputRoot '_fractions'];
fractions       =   envi2matlab(fractionsFile);
% % Flip to account for geo-referencing
fractions       =   flipud(fractions);
% Reset nodata and unmodeled pixels to NaN
fractions(fractions==0)=   NaN;
% % Apply mask
% fractions       =   applyMask(fractions,waterMask);


%% Display fraction images for the two taxa in the mixture, along with water
figure
nFrac   =   size(fractions,3);
titles  =   [Taxa(iTaxa2mix1); Taxa(iTaxa2mix2); 'Water'];
tiledlayout(1,nFrac);
for i =1:nFrac
    nexttile
    imagesc(fractions(:,:,i),'XData',R.XWorldLimits,'YData',R.YWorldLimits);
    axis image; axis xy; axis off;
    set(gca,'clim',[0 1]);
    crameri('-bamako');
    if i == 2
        c               =   colorbar('southOutside');
        c.Label.String  =   'End member fraction';
    end
    title(titles{i});
end


%% Compute differences between known mixture and MESMA output
inputMix    =   flipud(cat(3,fracTaxa1img,fracTaxa2img,fracWaterImg));
fracDiff    =   inputMix - fractions;
% figure
% imshow(fracDiff(:,:,1),[]);
% colorbar


%% Display fractions as a color composite: taxa1 as red, taxa2 as green, water as blue
figure
tiledlayout(1,3)
nexttile
imshow(inputMix,[0 1])
title({'Simulated mixture:',['R:' Taxa{iTaxa2mix1} ' G:' Taxa{iTaxa2mix2} ...
       ' B:Water']})
axis xy
nexttile
imshow(fractions,[0 1])
title({'MESMA output:',['R:' Taxa{iTaxa2mix1} ' G:' Taxa{iTaxa2mix2} ...
       ' B:Water']})
axis xy
nexttile
imshow(fracDiff,[0 1])
title({'Fraction differences:',['R:' Taxa{iTaxa2mix1} ' G:' Taxa{iTaxa2mix2} ...
       ' B:Water']})


%% *MESMA_rmse file contains the model's RMSE
% The model’s RMSE [nb of bands = 1]
% Value of unmodeled pixels in output: 9999
% Value of pixels with no data in output: 9998
rmseFile   =   [outputRoot '_rmse'];
rmse       =   envi2matlab(rmseFile);
% Flip to account for geo-referencing
rmse       =   flipud(rmse);
% Reset nodata pixels to NaN
rmse(rmse==9998)=   NaN;
% Reset unmodeled pixels to NaN;
rmse(rmse==9999)=   NaN;


%% Display RMSE image
figure
imagesc(rmse,'XData',R.XWorldLimits,'YData',R.YWorldLimits,...
        'alphadata',~isnan(rmse));
axis image; axis xy; axis off; 
colormap(gray)
colorbar
title('MESMA RMSE')


%% Summarize fraction errors and RMSE values
% Compute the mean and standard deviation of the fraction errors for each layer
% (i.e., end member) over all of the mixtures (i.e., the two horizontal
% dimensions of the fraction difference image)
meanFracError   =   squeeze(mean(fracDiff,[1 2]));
    % These are on the order of 1E-9, so tiny, essentially rounding error
for i = 1:nFrac
    tmp             =   fracDiff(:,:,i);
    stdFracError(i) =   std(tmp(:)); %#ok<SAGROW>
end
clear i tmp

% Compute mean and standard deviation of RMSE 
meanRmse        =   mean(rmse(:));  % Also tiny: 1.0760e-09
stdRmse         =   std(rmse(:));   % Also tiny: 4.5410e-10


%% For this simple, ideal test case, MESMA exactly reproduces known input fractions


%% Save what we have so far
clear ans i t c
save HABsimSpec.mat



%% NOTES FOR NEXT TIME ...
% Work on formalizing this prototyping code into functions
% Also for import and display of results
% Quantify reproduction of input fractions
% Experiment with other combinations of taxa based on ns3 scores
% Revisit OWASCO now that we know we can use a shade end member to represent
% water by calling the python implementation of MESMA