function []=make_pvd_file(pathname,filename,np,nt,ext)

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

f = fullfile(pathname,filename);
filename = strcat(f,'.pvd');
fmesh = fopen(filename,'w',endian_matlab);
fprintf(fmesh,'<?xml version="1" ?>\n');

fprintf(fmesh,['<VTKFile type="Collection" version="0.1" byte_order="',endian_paraview,'">\n']);
fprintf(fmesh,'<Collection>\n');

for i=1:np
    for j=0:nt-1
        fprintf(fmesh,'<DataSet timestep="%u" group="" part="%u" ',j,i);
        fprintf(fmesh,'file="%s_%u_%u.',f,i,j);
        fprintf(fmesh,ext);
        fprintf(fmesh,'"/>\n');
    end
end

fprintf(fmesh,'</Collection>\n');
fprintf(fmesh,'</VTKFile>\n');

fclose(fmesh);
