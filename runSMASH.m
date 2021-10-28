function outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
                               classNameField,imgSubFile,shadeFile,...
                               imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
                               minShadeEmFrac,maxShadeEmFrac,maxRMSE)
% Specify parameters and perform MESMA of algal taxa as part of HAB remote sensing workflow
%
%% runSMASH.m:
%   Set up input parameters and perform Multiple End Member Spectral Mixture
%   Analysis (MESMA) of algal taxa as the core of a workflow for remote sensing
%   of Harmful Algal Blooms (HABs).  This function takes the outputs from
%   preSMASH.m and the user-specified MESMA parameters as inputs and
%   programmatically generates a *.bat file to run a Python implementation of
%   MESMA.  The Python routine is called from within MATLAB and the MESMA
%   outputs are written to the specified directory.  MESMA results can be
%   imported and processed using postSMASH.m.
%
%% SYNTAX:
%   outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
%                         classNameField,imgSubFile,shadeFile,...
%                         imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
%                         minShadeEmFrac,maxShadeEmFrac,maxRMSE);
%                            
%% INPUTS:
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
%   maxRMSE:        Maximum RMSE constraint used in MESMA model
%
%% OUTPUTS:
%   outputRoot:     Root file name for MESMA outputs, with a time stamp used to
%                   identify the MESMA run
%
%% NOTES:
% > Code contains comments with inputs used for Lake Owasco prototype
% > Takes output from preSMASH.m as input
% > Programatically calls the Python implementation of MESMA, which involves
%   some hard-wired paths and code
% > To import MESMA results, use postSMASH.m
% > For information on the core algorithm, see the Python MESMA documentation at
%   https://mesma.readthedocs.io/en/latest/userguide/mesma_cli.html
%
%% FUNCTION SUMMARY:
%   outputRoot = runSMASH(dataDir,siteCode,dateCode,libraryFile,...
%                         classNameField,imgSubFile,shadeFile,...
%                         imgScaleFactor,libScaleFactor,minEmFrac,maxEmFrac,...
%                         minShadeEmFrac,maxShadeEmFrac,maxRMSE);

%% CREDITS:
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 06/18/2021
% Also see prototyping code in: C:\Users\cjl\OneDrive - DOI\HABs\SMASH\SMASHowasco.m
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\runSMASH.m


%% EXAMPLE INPUTS: Set parameters we will pass to the MESMA core algorithm
% % Name of input library
% libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\AlgaeTaxaOnly.sli';
% % Metadata field with end-member names
% classNameField  =   'Class';
% % Name of image file
% imgSubFile      =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Owasco\OW20200823specSub.tif';
% % Scale factors for library and image
% libScaleFactor  =   1;
% imgScaleFactor  =   10000;
% % Constraints
% %   Minimum end member fraction
% minEmFrac       =   -0.05;
% %   Maximum end member fraction
% maxEmFrac       =   1.05;
% %   Minimum shade end member fraction
% minShadeEmFrac  =   -0.05;
% %   Maximum shade end member fraction
% maxShadeEmFrac  =   1.05;
% %   Maximum RMSE
% maxRMSE         =   0.025;
% % Name of shade end member file
% shadeFile       =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\DetroitLakeWaterOnly.sli';


%% Programmatically generate new magic.bat file with actual call to mesma embedded
% Set up the actual call to the mesma routine itself first so we can embed it in
% the .bat file we're writing here
siteDateCode=   [siteCode dateCode];
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


%% Now call the new .bat file to run the mesma
% Specify the command as a string and then call it using the system function
command =   [batFileName ' &'];
status  =   system(command);
if status == 0
    disp('MESMA executed successfully, proceed to postSMASH.m to import results ...')
else
    disp('WARNING: MESMA did not execute as intended ...')
end