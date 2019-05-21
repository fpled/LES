function [] = make_pvd_file(pathname,filename,np,nt,ext)
% function [] = make_pvd_file(pathname,filename,np,nt,ext)

if nargin<1 || isempty(pathname)
    pathname='.';
end
if nargin<2 || isempty(filename)
    filename='paraview';
end
if nargin<3 || isempty(np)
    np=1;
end
if nargin<4 || isempty(nt)
    nt=1;
end
if nargin<5 || isempty(ext)
    ext='vtu';
end

[~,~,endian]=computer;
if endian=='L'
    endian_matlab='ieee-le';
    endian_paraview='LittleEndian';
else
    endian_matlab='ieee-be';
    endian_paraview='BigEndian';
end

f = filename;
filename = fullfile(pathname,strcat(filename,'.pvd'));
fid = fopen(filename,'w',endian_matlab);
fprintf(fid,'<?xml version="1.0" ?>\n');

fprintf(fid,['<VTKFile type="Collection" version="0.1" byte_order="',endian_paraview,'">\n']);
fprintf(fid,'\t <Collection>\n');

for i=1:np
    for j=0:nt-1
        fprintf(fid,'\t\t <DataSet timestep="%u" group="" part="%u" ',j,i);
        fprintf(fid,'file="%s_%u_%u.',f,i,j);
        fprintf(fid,ext);
        fprintf(fid,'"/>\n');
    end
end

fprintf(fid,'\t </Collection>\n');
fprintf(fid,'</VTKFile>\n');

fclose(fid);
