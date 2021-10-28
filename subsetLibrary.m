function librarySub = subsetLibrary(libraryFile)
% Select taxa to include in a subset of the algal spectral library
% 
%% subsetLibrary.m:
%   Helper function for selecting algal taxa to include in a subset of the
%   spectral library. Intended for use within the SMASH framework. Note that
%   actual spectral library data file is hard-wired into code. The user is
%   prompted choose select the taxa ID numbers to be included in the subset of
%   the spectral library based on a plot of the spectra and a legend that lists
%   the ID numbers for each taxa.
% 
%% SYNTAX: 
%   librarySub  =   subsetLibrary(libraryFile);
% 
%% INPUTS:
%   libraryFile:    String specifying name of the full algal spectral library to
%                   be imported, displayed, and subset. Note that actual
%                   spectral library data file is hard-wired into code.
%
%% OUTPUTS:
%   librarySub:     Subset of the original spectral library containing only the
%                   selected taxa; the corresponding metadata *.CSV file is also
%                   created. The library file created is named LibrarySubset.csv
%                   and is saved in the current working directory.
% 
%% NOTES:
%   Developed for use with the Upper Klamath Lake data set to exclude one taxa
%   and repeat the smash analysis. See also preSmash, runSmash, and postSmash
%   functions. Note that actual spectral library data file is hard-wired into
%   code.
%
%% FUNCTION SUMMARY:
%   librarySub  =   subsetLibrary(libraryFile); 

%% CREDITS:
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 07/13/2021, 7/29/2021
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\subsetLibrary.m


%% First import and display the full library
% cd('C:\Users\cjl\OneDrive - DOI\HABs\SMASH')
% Import the library
% libraryFile     =   'C:\Users\cjl\OneDrive - DOI\HABs\SMASH\Library\Algae+Water.sli';
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


%% Actual spectral library data file hard-wired into code
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


%% Subset the library to include only the selected algal taxa
% Prompt user to specify indices of taxa to be included in subset library
iTaxa       =   input('Enter a vector of spectrum numbers of taxa to include in subset library: ');
librarySub  =   library(iTaxa,:);
% Best way to proceed will be to just write a new ENVI header using the code in
% the new function specLibHeader.m after calling matlab2envi to actually write
% out the data
matlab2envi(librarySub,'LibrarySubset.sli',[],[],[],wvlSub);
Header      =   specLibHeader(librarySub,'LibrarySubset.sli',Taxa(iTaxa),wvlSub); %#ok<*NASGU>
% Also make the corresponding csv metadata file
Metadata    =   table(Taxa(iTaxa),...
                      Taxa(iTaxa),...
                      'VariableNames',{'Name','Class'});
writetable(Metadata,'LibrarySubset.csv')