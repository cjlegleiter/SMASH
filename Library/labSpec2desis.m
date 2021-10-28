%% SMASH SPECTRAL LIBRARY RESAMPLING: LAB SPECTRA TO DESIS BANDS
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 05/08/2021
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\labSpec2desis.m


%% Info on lab spectra
%   Tyler King provided lab spectra: 
%       C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Hab_Spectra_Corrected.zip
% 	Each *.csv file is named for the algal taxa and has two columns of data
% 	starting in the fifth row: wavelength in nm and "corrected reflectance"
%   The first row of each file has the taxa name


%% Unzip lab spectra and read each csv file into a MATLAB structure array
cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
unzip('C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Hab_Spectra_Corrected.zip',...
      'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\labSpectra')
fileList    =   dir('labSpectra\*.csv');
LabSpectra  =   struct('taxa',[],'wvl',[],'rfl',[]);

% Used the GUI import tool for one of the files and asked it to write a script,
% which provided the following code
    % Set up the Import Options and import the data
    opts    =   delimitedTextImportOptions("NumVariables", 2, "Encoding", "UTF-8");
    % Specify range and delimiter
    opts.DataLines          =   [5, Inf];
    opts.Delimiter          =   ",";
    % Specify column names and types
    opts.VariableNames      =   ["wvl","rfl"];
    opts.VariableTypes      =   ["double","double"];
    % Specify file level properties
    opts.ExtraColumnsRule   =   "ignore";
    opts.EmptyLineRule      =   "read";

% Now loop over the files and read in each one in turn    
for i = 1:length(fileList)
    % Get the name
    name    =   extractBefore(fileList(i).name,'.csv');
    % Import the data
    tmp     =   readtable(fullfile(fileList(i).folder,fileList(i).name),opts);
    % Allocate to structure array
    LabSpectra(i).taxa  =   name;
    LabSpectra(i).wvl   =   tmp.wvl;
    LabSpectra(i).rfl   =   tmp.rfl;
end


%% FWHM info for hyperspectral microscope from Terry Slonecker
% The guys at Surface Optics Corp that made the SOC 710 camera on David Allen’s
% hyperspectral microscope say the FWHM for that instrument is 4.69 nm.
fwhm    =   4.69;
% Assume this is the same for all bands
for i = 1:length(LabSpectra)
    LabSpectra(i).fwhm  =   ones(size(LabSpectra(i).wvl))*fwhm;
end
% Reorder fields
LabSpectra  =   orderfields(LabSpectra,{'taxa','wvl','fwhm','rfl'});


%% Write spectra out to a text file for input to ENVI, including FWHM column
% Conver structure array to table
LabSpecTable    =   struct2table(LabSpectra);
% Get the list of taxa as a cell array of characters; note transpose to make
% this into a row rather than a column
taxaList        =   LabSpecTable.taxa';
% Set up an array to hold list of wavelengths, fwhm, and then reflectance values
% for each taxa
labSpecData     =   zeros(length(LabSpectra(1).wvl),length(taxaList)+2);
% Make the wavelengths the first column
labSpecData(:,1)=   LabSpectra(1).wvl;
% Make the fwhm the second column
labSpecData(:,2)=   LabSpectra(1).fwhm;
% Loop over the taxa to populate the remaining columns
for i = 1:length(taxaList)
    % Also apply scale factor of 10000 to convert to reflectance between 0 and 1
    labSpecData(:,i+2)  =   LabSpectra(i).rfl/10000;
end
% Plot to confirm
figure
plot(labSpecData(:,1),labSpecData(:,3:end))
xlabel('Wavelength (nm)')
ylabel('Reflectance')
legend(taxaList)
title('HAB Spectral Library')
% Make header row with wavelength in first column, then list of taxa
headerRow   =   ["Wavelength" taxaList];
% Write out the data to a text file, starting with the header
outFile     =   'HABlabSpec.csv';
writematrix(headerRow,outFile)
% Now append the actual data to the same file
writematrix(labSpecData,outFile,'WriteMode','append')
% Also make a separate text file with just the names, which is what ENVI expects
outFile     =   'HABlabSpecNames.txt';
writecell(taxaList',outFile)


%% Next steps implemented in ENVI
% 	• Import the lab spectra to ENVI
% 		○ Toolbox > Spectral Libraries > Spectral Library Builder
% 		○ Input Spectral Wavelength From: ASCII File
% 		○ Select HABlabSpec.csv
% 		○ Wavelength Column = 1, FWHM column = 2, Wavelength Units =
% 		Nanometers, Y Scale Factor = 1
% 		○ Spectral Library Builder window opens
% 		○ Import > From ASCII File; select HABlabSpec.csv
% 		○ X Axis Column = 1; Select Y Axis Columns: 3-15; Wavelength Units =
% 		Nanometers; Y Scale Factor = 1
% 		○ Options > Import spectrum names from ASCII; select HABlabSpecNames.txt
% 		○ Click Select All and then Plot to plot the spectra, then turn on the
% 		legend under Options in the plot window
% 		○ In the Spectral Library Builder window, choose File>Save spectra
% 		as>Spectral library file and populate the dialog as follows, then save
% 		as HABlabSpec.sli









%% Set up DESIS spectral response function
% Set up filter function for DESIS bands using the spectral response functions
% for each band that Tyler provided: C:\Users\cjl\OneDrive -
% DOI\HABs\SMASH\desis_bands.zip


%% Unzip desis band info folder and import response function for each band
inFile      =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\desis_bands.zip';
unzip(inFile,'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\desis_bands')
fileList    =   dir('.\desis_bands\desis_bands\desis_response_band_*.csv');
% Set up structure array to store bands
DesisBands  =   struct('band',[],'wvl',[],'response',[]);
% Used the GUI import tool for one of the files and asked it to write a script,
% which provided the following code
    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 3);
    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";
    % Specify column names and types
    opts.VariableNames = ["Var1", "wavelength", "response"];
    opts.SelectedVariableNames = ["wavelength", "response"];
    opts.VariableTypes = ["char", "double", "double"];
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    % Specify variable properties
    opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");
    
    
%% Get list of DESIS band centers Tyler extracted from XML
inFile          =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\desis_bands\desis_bands\desis_bands_from_xml.csv';
desisBandCenters=   readtable(inFile);
desisBandCenters=   desisBandCenters.central_wl;


%% Loop over the bands, each of which has its own file
for i = 1:length(fileList)
    % Import the data
    fileName                =   fullfile(fileList(i).folder,fileList(i).name);
    tmp                     =   readtable(fileName,opts);
    % Be careful with file names and band numbers and how they are sorted in the
    % fileList structure array: for bands less than 10 the d_ precedes the band
    % number, for bands between 10 and 99 the _ precedes the band number, and
    % for bands 100 and above we can just use the band number
    tmpStr                  =   fileName(end-6:end-4);
    if strcmp(tmpStr(1),'_')
        tmpStr  =   tmpStr(2:end);
    elseif strcmp(tmpStr(1),'d')
        tmpStr  =   tmpStr(end);
    end
    DesisBands(i).band      =   str2double(tmpStr);
    DesisBands(i).wvl       =   tmp.wavelength;
    DesisBands(i).response  =   tmp.response;
end
% Sort the structure array to make sure the bands are increasing from 1 to 235
[~,iSort]   =   sort([DesisBands.band]);
DesisBands  =   DesisBands(iSort);
figure
hold on
for i = 1:length(DesisBands)
    plot(DesisBands(i).wvl,DesisBands(i).response);
    text(DesisBands(i).wvl(25),DesisBands(i).response(25),num2str(i));
    drawnow
    pause(0.1)    
end


%% Loop over bands and get response function wavelengths for all bands
% First we need the wavelengths in each band's response function
allWvl      =   [];
for i = 1:length(DesisBands)
    % Note the NaN at the beginning of each list of wavelengths, which we should
    % just get rid of now
    wvlTmp  =   DesisBands(i).wvl(2:end);
    allWvl  =   [allWvl; wvlTmp]; %#ok<AGROW>
end
% Sort so the wavelengths are monotonically increasing
allWvl      =   sort(allWvl);


%% Now to populate the response function matrix, with wavelengths as first column
allResp     =   zeros(length(allWvl),length(DesisBands)+1);
allResp(:,1)=   allWvl;
for i = 1:length(DesisBands)
    % Remove NaN at the beginning
    wvlTmp          =   DesisBands(i).wvl(2:end);
    respTmp         =   DesisBands(i).response(2:end);
    % Also normalize so values range from 0 to 1
    respTmp         =   respTmp/max(respTmp);
    % Interpolate response function values for this band at the wavelengths for
    % all of the response function wavelengths from all of the bands
    respInt         =   interp1(wvlTmp,respTmp,allWvl,'linear',0);
    % Populate array with interpolated values of response function for this
    % band, but don't forget the wavelengths themselves are the first column
    allResp(:,i+1)  =   respInt;
end
figure
plot(allWvl,allResp(:,2:end))
xlabel('Wavelength (nm)')
ylabel('Normalized spectral response')
title('DESIS band filter functions')


%% Discard bands beyond range of wavelengths in lab spectra
labSpecWvl      =   LabSpecTable.wvl{1};
maxWvl          =   max(labSpecWvl)+mean(diff(labSpecWvl))/2;
iKeep           =   find(desisBandCenters<maxWvl);
% Don't forget that the first column is the wavelengths
allRespSub      =   allResp(:,[1; iKeep+1]);


%% Output as a text file we can read in to ENVI as a spectral library
outFile     =   'DesisBandResponse.csv';
writematrix(allRespSub,outFile)


%% Also make a separate text file with just the names, which is what ENVI expects
outFile     =   'DesisBandNames.txt';
BandNames   =   cell(size(DesisBands(iKeep)));
for i = 1:length(DesisBands(iKeep))
   BandNames{i}     =   ['DesisBand' num2str(i)];
end
writecell(BandNames',outFile)


%% Next steps in ENVI; see OneNote HABs+WQ / SPECTRAL LIBRARY RESAMPLING
% 	• Import the DESIS response function to ENVI
% 		○ Toolbox > Spectral Libraries > Spectral Library Builder
% 		○ Input Spectral Wavelength From: First input spectrum
% 		○ In Spectral Library Builder window, select Import>from ASCII file
% 		○ Select DesisBandResponse.csv
% 		○ X axis column = 1, Select all remaining columns (2:159) for Y Axis, Wavelength Units = Nanometers, Y Scale Factor = 1
% 		○ Options > Import spectrum names from ASCII; select DesisBandNames.txt
% 		○ Click Select All and then File>Save spectra as>Spectral library file
% 		
% 		
% 		
% 		Screen clipping taken: 5/9/2021 3:32 PM
% 		
% 		○ Save as C:\Users\cjl\OneDrive - DOI\HABs\SMASH\DesisBandResponse.sli
% 		○ The spectral library viewer will open and you can plot the response functions for specific bands
% 	• To perform the actual spectral resampling …
% 		○ From the ENVI toolbox, select Spectral>Spectral Libraries>Spectral Library Resampling
% 		○ Select HABlabSpec.sli as the input file
% 		○ Next, Resample Wavelength to: User Defined Filter Function, output to C:\Users\cjl\OneDrive - DOI\HABs\SMASH\HABlab2desis.sli
% 		○ For the Input Filter Function Spectral Library, choose HABlabSpec.sli
% 		○ You might get a warning about some of the points in the output wavelength being outside the range of the input wavelength and not being interpolated, but click OK
% 		○ You can plot the original and convolved spectra with the Spectral Library Viewer
% 		○ 
% 		○ To export the convolved spectra to a text file, click on each one in turn to add them to the Spectral Library Viewer plot window and then click on the Export drop down and select ASCII
% 		○ Saved as C:\Users\cjl\OneDrive - DOI\HABs\SMASH\HABlab2desis.txt
% 		○ There's a header with plotting info, but you should be able to import this into other software


%% Save all this to a .mat file
save labSpec2desis.mat