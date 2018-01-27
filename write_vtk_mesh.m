function []=write_vtk_mesh(M,u,phase,pathname,filename,part,time,binary_output)

if nargin<6
    part=1;
end
if nargin<7
    time=0;
end
if nargin<8
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

filename=fullfile(pathname,strcat(filename,'_',num2str(part),'_',num2str(time),'.vtu'));
fmesh = fopen(filename,'w',endian_matlab);
fprintf(fmesh,'<?xml version="1" ?>\n');

fprintf(fmesh,['<VTKFile type="UnstructuredGrid" version="0.1" byte_order="',endian_paraview,'">\n']);
fprintf(fmesh,'\t <UnstructuredGrid>\n');
fprintf(fmesh,'\t\t <Piece NumberOfPoints="%u" NumberOfCells="%u">\n',nnode,nelem);

if binary_output
    
    % POINT DATA
    fprintf(fmesh,'\t\t\t <PointData scalars="scalar"> \n');
    fprintf(fmesh,'\t\t\t\t <DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(u);
    fprintf(fmesh,'\t\t\t\t <DataArray type="Float32" Name="phase" NumberOfComponents="1" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(phase);
    fprintf(fmesh,'\t\t\t </PointData> \n');
    
    % CELL DATA
    fprintf(fmesh,'\t\t\t <CellData> \n');
    fprintf(fmesh,'\t\t\t </CellData> \n');
    
    % POINTS
    fprintf(fmesh,'\t\t\t <Points>\n');
    fprintf(fmesh,'\t\t\t\t <DataArray type="Float32" NumberOfComponents="3" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(node);
    fprintf(fmesh,'\t\t\t </Points>\n');
    
    % CELLS
    fprintf(fmesh,'\t\t\t <Cells>\n');
    fprintf(fmesh,'\t\t\t\t <DataArray type="Int32" Name="connectivity" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(connectivity);
    fprintf(fmesh,'\t\t\t\t <DataArray type="Int32" Name="offsets" format="appended" offset="%u" />\n',offset);
    offset=offset+4+4*numel(offsets);
    fprintf(fmesh,'\t\t\t\t <DataArray type="Int8" Name="types" format="appended" offset="%u" />\n',offset);
    fprintf(fmesh,'\t\t\t </Cells>\n');
    
    % APPENDED DATA
    fprintf(fmesh,'\t\t </Piece>\n');
    fprintf(fmesh,'\t </UnstructuredGrid> \n');
    fprintf(fmesh,'\t <AppendedData encoding="raw"> \n _');
    
    % DISPLACEMENT U
    fwrite(fmesh,4*numel(u),'uint32');
    fwrite(fmesh,u,'float32');
    
    % PHASE
    fwrite(fmesh,4*numel(phase),'uint32');
    fwrite(fmesh,phase,'float32');
    
    % NODES
    fwrite(fmesh,4*numel(node),'uint32');
    fwrite(fmesh,node','float32');
    
    % ELEMS
    fwrite(fmesh,4*numel(connectivity),'uint32');
    fwrite(fmesh,connectivity'-1,'int32');
    fwrite(fmesh,4*numel(offsets),'uint32');
    fwrite(fmesh,offsets,'int32');
    fwrite(fmesh,numel(types),'uint32');
    fwrite(fmesh,types,'int8');
    
    fprintf(fmesh,'\n%s\n','</AppendedData>');
else
    
    % POINT DATA
    fprintf(fmesh,'\t\t\t <PointData scalars="scalar"> \n');
    fprintf(fmesh,'\t\t\t\t <DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="ascii">\n');
    fprintf(fmesh,'%e \n',u);
    fprintf(fmesh,'\t\t\t\t </DataArray>\n');
    fprintf(fmesh,'\t\t\t\t <DataArray type="Float32" Name="phase" NumberOfComponents="1" format="ascii">\n');
    fprintf(fmesh,'%e \n',phase);
    fprintf(fmesh,'\t\t\t\t </DataArray>\n');
    fprintf(fmesh,'\t\t\t </PointData> \n');
    
    % CELL DATA
    fprintf(fmesh,'\t\t\t <CellData> \n');
    fprintf(fmesh,'\t\t\t </CellData> \n');
    
    % POINTS
    fprintf(fmesh,'\t\t\t <Points>\n');
    fprintf(fmesh,'\t\t\t\t <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
    fprintf(fmesh,'%f \n',node');
    fprintf(fmesh,'\t\t\t\t </DataArray>\n');
    fprintf(fmesh,'\t\t\t </Points>\n');
    
    % CELLS
    fprintf(fmesh,'\t\t\t <Cells>\n');
    fprintf(fmesh,'\t\t\t\t <DataArray type="Int32" Name="connectivity" format="ascii">\n');
    fprintf(fmesh,'%u \n',connectivity'-1);
    fprintf(fmesh,'\t\t\t\t </DataArray>\n');
    fprintf(fmesh,'\t\t\t\t <DataArray type="Int32" Name="offsets" format="ascii">\n');
    fprintf(fmesh,'%u \n',offsets);
    fprintf(fmesh,'\t\t\t\t </DataArray>\n');
    fprintf(fmesh,'\t\t\t\t <DataArray type="Int8" Name="types" format="ascii">\n');
    fprintf(fmesh,'%u \n',types);
    fprintf(fmesh,'\t\t\t\t </DataArray>\n');
    fprintf(fmesh,'\t\t\t </Cells>\n');
    
    % END VTK FILE
    fprintf(fmesh,'\t\t </Piece>\n');
    fprintf(fmesh,'\t </UnstructuredGrid> \n');
end
fprintf(fmesh,'</VTKFile> \n');

fclose(fmesh);

end
