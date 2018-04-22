function []=write_vtk_mesh(M,u,C,tauConv,tauDiff,tauSurf,tauInterf,pathname,filename,part,time,binary_output)

if nargin<8 || isempty(pathname)
    pathname='.';
end
if nargin<9 || isempty(filename)
    filename='paraview';
end
if nargin<10 || isempty(part)
    part=1;
end
if nargin<11 || isempty(time)
    time=0;
end
if nargin<12 || isempty(binary_output)
    binary_output=1;
end

nelem=M.nbelem;
nnode=M.nbnode;
connectivity=[];
offsets=[];
types=[];
for i=1:M.nbgroupelem
    elem = M.groupelem{i};
    nbnode = getnbnode(elem);
    nbelem = getnbelem(elem);
    co = getconnec(elem)';
    [~,co] = ismember(co,getnumber(M.node));
    connectivity = [connectivity ; co(:)];
    if isempty(offsets)
        coeffs=0;
    else
        coeffs=offsets(size(offsets,2));
    end
    offsets = [offsets coeffs+nbnode*(1:nbelem)];
    type = vtkelemtype(elem);
    types=[types repmat(type,1,getnbelem(elem))];
end

x = getcoord(M.node)';
node = x(:);

[~,~,endian]=computer;
if endian=='L'
    endian_matlab='ieee-le';
    endian_paraview='LittleEndian';
else
    endian_matlab='ieee-be';
    endian_paraview='BigEndian';
end
offset=0;

filename = fullfile(pathname,strcat(filename,'_',num2str(part),'_',num2str(time),'.vtu'));
fid = fopen(filename,'w',endian_matlab);
fprintf(fid,'<?xml version="1.0" ?>\n');

fprintf(fid,['<VTKFile type="UnstructuredGrid" version="0.1" byte_order="',endian_paraview,'">\n']);
fprintf(fid,'\t <UnstructuredGrid>\n');
fprintf(fid,'\t\t <Piece NumberOfPoints="%u" NumberOfCells="%u">\n',nnode,nelem);

if binary_output
    
    % POINT DATA
    fprintf(fid,'\t\t\t <PointData scalars="scalar"> \n');
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(u);
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="phase" NumberOfComponents="1" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(C);
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="tau conv" NumberOfComponents="3" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(tauConv);
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="tau diff" NumberOfComponents="3" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(tauDiff);
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="tau surf" NumberOfComponents="3" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(tauSurf);
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="tau interf" NumberOfComponents="1" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(tauInterf);
    fprintf(fid,'\t\t\t </PointData> \n');
    
    % CELL DATA
    fprintf(fid,'\t\t\t <CellData> \n');
    fprintf(fid,'\t\t\t </CellData> \n');
    
    % POINTS
    fprintf(fid,'\t\t\t <Points>\n');
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" NumberOfComponents="3" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(node);
    fprintf(fid,'\t\t\t </Points>\n');
    
    % CELLS
    fprintf(fid,'\t\t\t <Cells>\n');
    fprintf(fid,'\t\t\t\t <DataArray type="Int32" Name="connectivity" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(connectivity);
    fprintf(fid,'\t\t\t\t <DataArray type="Int32" Name="offsets" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(offsets);
    fprintf(fid,'\t\t\t\t <DataArray type="Int8" Name="types" format="appended" offset="%u" />\n',offset);
    fprintf(fid,'\t\t\t </Cells>\n');
    
    % APPENDED DATA
    fprintf(fid,'\t\t </Piece>\n');
    fprintf(fid,'\t </UnstructuredGrid> \n');
    fprintf(fid,'\t <AppendedData encoding="raw"> \n _');
    
    % DISPLACEMENT U
    fwrite(fid,4*numel(u),'uint32');
    fwrite(fid,u,'float32');
    
    % PHASE C
    fwrite(fid,4*numel(C),'uint32');
    fwrite(fid,C,'float32');
    
    % CONVECTION STRESS TAU CONV
    fwrite(fid,4*numel(tauConv),'uint32');
    fwrite(fid,tauConv,'float32');
    
    % DIFFUSION (VISCOSITY) STRESS TAU DIFF
    fwrite(fid,4*numel(tauDiff),'uint32');
    fwrite(fid,tauDiff,'float32');
    
    % SURFACE TENSION (CAPILLARY) STRESS TAU SURF
    fwrite(fid,4*numel(tauSurf),'uint32');
    fwrite(fid,tauSurf,'float32');
    
    % INTERFACE TRACKING STRESS TAU INTERF
    fwrite(fid,4*numel(tauInterf),'uint32');
    fwrite(fid,tauInterf,'float32');
    
    % NODES
    fwrite(fid,4*numel(node),'uint32');
    fwrite(fid,node','float32');
    
    % ELEMS
    fwrite(fid,4*numel(connectivity),'uint32');
    fwrite(fid,connectivity'-1,'int32');
    fwrite(fid,4*numel(offsets),'uint32');
    fwrite(fid,offsets,'int32');
    fwrite(fid,numel(types),'uint32');
    fwrite(fid,types,'int8');
    
    fprintf(fid,'\n%s\n','</AppendedData>');
else
    
    % POINT DATA
    fprintf(fid,'\t\t\t <PointData scalars="scalar"> \n');
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="ascii">\n');
    fprintf(fid,'%e \n',u);
    fprintf(fid,'\t\t\t\t </DataArray>\n');
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="phase" NumberOfComponents="1" format="ascii">\n');
    fprintf(fid,'%e \n',C);
    fprintf(fid,'\t\t\t\t </DataArray>\n');
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="tau conv" NumberOfComponents="3" format="ascii">\n');
    fprintf(fid,'%e \n',tauConv);
    fprintf(fid,'\t\t\t\t </DataArray>\n');
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="tau diff" NumberOfComponents="3" format="ascii">\n');
    fprintf(fid,'%e \n',tauDiff);
    fprintf(fid,'\t\t\t\t </DataArray>\n');
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="tau surf" NumberOfComponents="3" format="ascii">\n');
    fprintf(fid,'%e \n',tauSurf);
    fprintf(fid,'\t\t\t\t </DataArray>\n');
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="tau Interf" NumberOfComponents="1" format="ascii">\n');
    fprintf(fid,'%e \n',tauInterf);
    fprintf(fid,'\t\t\t\t </DataArray>\n');
    fprintf(fid,'\t\t\t </PointData> \n');
    
    % CELL DATA
    fprintf(fid,'\t\t\t <CellData> \n');
    fprintf(fid,'\t\t\t </CellData> \n');
    
    % POINTS
    fprintf(fid,'\t\t\t <Points>\n');
    fprintf(fid,'\t\t\t\t <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
    fprintf(fid,'%f \n',node');
    fprintf(fid,'\t\t\t\t </DataArray>\n');
    fprintf(fid,'\t\t\t </Points>\n');
    
    % CELLS
    fprintf(fid,'\t\t\t <Cells>\n');
    fprintf(fid,'\t\t\t\t <DataArray type="Int32" Name="connectivity" format="ascii">\n');
    fprintf(fid,'%u \n',connectivity'-1);
    fprintf(fid,'\t\t\t\t </DataArray>\n');
    fprintf(fid,'\t\t\t\t <DataArray type="Int32" Name="offsets" format="ascii">\n');
    fprintf(fid,'%u \n',offsets);
    fprintf(fid,'\t\t\t\t </DataArray>\n');
    fprintf(fid,'\t\t\t\t <DataArray type="Int8" Name="types" format="ascii">\n');
    fprintf(fid,'%u \n',types);
    fprintf(fid,'\t\t\t\t </DataArray>\n');
    fprintf(fid,'\t\t\t </Cells>\n');
    
    % END VTK FILE
    fprintf(fid,'\t\t </Piece>\n');
    fprintf(fid,'\t </UnstructuredGrid> \n');
end
fprintf(fid,'</VTKFile> \n');

fclose(fid);

end
