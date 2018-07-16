function fout = read3b42(filein)

%fout = read3b42rt(filein)
%
%Function to read the 3B42RT binary data files obtained from
%    ftp://trmmopen.gsfc.nasa.gov/pub/merged/mergeIRMicro
%    and write the data out to a MATLAB struct.
%
%
%INPUT:
%filein: the name of the HDF file to read
%
%OUTPUT:
%fout: the output struct
%
%------------------------------------------
%Brandon Kerns, RSMAS
%bkerns@orca.rsmas.miami.edu
%3.24.2009
%------------------------------------------
%CHANGES:
%
%
%------------------------------------------

if ischar(filein)


fid = fopen(filein,'r','b')  ;


fout.lon = -179.875:0.25:179.875 ;
fout.lat = -49.875:0.25:49.875 ;

fout.precip = double(hdfread(filein,'precipitation'))' ;
fout.precip_ir = double(hdfread(filein,'IRprecipitation'))' ;
fout.precip_mw = double(hdfread(filein,'HQprecipitation'))' ;
fout.errors = double(hdfread(filein,'relativeError'))' ;
fout.source = double(hdfread(filein,'satPrecipitationSource'))' ;
fout.source_description = char({...
    '0   no observation   ',...
    '1   AMSU   ',...
    '2   TMI   ',...
    '3   AMSR   ',...
    '4   SSMI   ',...
    '5   SSMI/S   ',...
    '6   MHS   ',...
    '7   TCI   ',...
    '30  AMSU/MHS average   ',...
    '31  Conical scanner average (includes TMI, AMSR, SSMI, and SSMI/S.)   ',...
    '50  IR   ',...
    'Add 100 to above if sampling is less than or equal to two pixels'});
fout.toffset = double(hdfread(filein,'satObservationTime'))' ;
fout.toffset_description = 'Minutes Relative to Nominal Time' ;

fclose(fid);

else

    'Input must be a string.'
    fout=''
	
end


