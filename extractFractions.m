function fracAtField = extractFractions(fieldLat,fieldLon,siteCode,dateCode,...
                                        epsgCode,fractions,Taxa,R)
% Extract end member fractions from MESMA output pixel at location of field sample
% 
%% extractFractions.m:
%   Post-processing function for extracting end member fractions from the MESMA
%   output image for the pixel at a specified location where a field sample was
%   collected.  These fraction estimates can then be compared to observed
%   biovolumes as a form of accuracy assessment.

%% SYNTAX:
% fracAtField   =   extractFractions(fieldLat,fieldLon,siteCode,dateCode,...
%                                    epsgCode,fractions,Taxa,R);
% 
%% INPUTS:
%   fieldLat:   latitude of the field sample location
%   fieldLon:   longitude of the field sample location
%   siteCode:   String with two-letter site code
%   dateCode:   String with 8-digit date code (YYYYMMDD)
%   epsgCode:   EPSG projection code for the image, used to convert the field
%               sampling location from lat, long to UTM
%   fractions:  Image with each band representing the MESMA fractions for each
%               of the algal taxa, plus an additional (last) band for the water
%               end member used as shade
%   Taxa:       Cell array of strings with names of algal taxa in the library 
%               that served as end members for MESMA
%   R:          Geo-referencing object for the fraction image
% 
%% OUTPUT:
%   fracAtField:(nTaxa+1) X 1 vector with end member fractions for each of the
%               nTaxa taxa in the library, plus the water fraction, extracted
%               from the MESMA output pixel at the field sample location
% 
%% FUNTION SUMMARY:
% fracAtField   =   extractFractions(fieldLat,fieldLon,siteCode,dateCode,...
%                                    epsgCode,fractions,Taxa,R);


%% CREDITS:
% Dr. Carl J. Legleiter, cjl@usgs.gov            
% Geomorphology and Sediment Transport Laboratory
% United States Geological Survey                
% 08/12/2021
% C:\Users\cjl\OneDrive - DOI\HABs\SMASH\extractFractions.m


%% Example inputs for Rattlesnake Point in UKL
% fieldLat        =   42.344214;
% fieldLon        =   -121.858296;


%% Project from lat, long to UTM 
[fieldX,fieldY] =   projfwd(projcrs(epsgCode),fieldLat,fieldLon);


%% Extract fractions for all taxa at this location
% fracAtField     =   zeros(length(Taxa)+1,1);
fracAtField     =   zeros(size(fractions,3));
for i = 1:length(fracAtField)
    tmp             =   impixel(R.XWorldLimits,R.YWorldLimits,...
                                fractions(:,:,i),fieldX,fieldY);
    fracAtField(i)  =   tmp(1);
end


%% Plot bar graph with fractions for each taxa + water
siteDateCode    =   [siteCode dateCode];
figure
bar(fracAtField)
ylabel('End member fraction')
% set(gca,'xticklabels',[Taxa; {'Water'}])
% Account for inclusion of unclassified in list of taxa
set(gca,'xticklabels',[Taxa(2:end); {'Water'}])
title([siteDateCode ': end member fractions for algal taxa'])