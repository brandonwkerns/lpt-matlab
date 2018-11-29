function fdata = read3b42rt_struct(filein)

%
%fdata = read3b42rt_struct(filein)
%
%Function to read the 3B42RT binary data files obtained from
%    ftp://trmmopen.gsfc.nasa.gov/pub/merged/mergeIRMicro
%    and return a matlab data struct.
%
%
%INPUT:
%filein: the name of the binary file to read
%fileout: the name of the HDF file to write to, default to filein.HDF
%
%OUTPUT:
%fdata, the matlab data struct.
%
%------------------------------------------
%Brandon Kerns, RSMAS
%bkerns@orca.rsmas.miami.edu
%8.12.2008
%------------------------------------------
%CHANGES:
%
%
%------------------------------------------

if ischar(filein)

if nargin < 2
	fileout = [filein,'.HDF'] ;
end

fid = fopen(filein,'r','b')  ;

%%%read header
header = fgets(fid, 2880) ;

fdata.precip = flipdim(reshape(fread(fid,691200,'int16'),1440,480)/100.0,2)' ;
fdata.errors = flipdim(reshape(fread(fid,691200,'int16'),1440,480)/100.0,2)' ;
fdata.source = flipdim(reshape(fread(fid,691200,'int8'),1440,480),2)' ;

fdata.lon = 0.125:0.25:359.875 ;
fdata.lat = fliplr(59.875:-0.25:-59.875) ;


fclose(fid) ;

else

	'Input must be a string.'
	fileout=''
	
end


