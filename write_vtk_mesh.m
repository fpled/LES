function [] = write_vtk_mesh(M,nodalfields,elemfields,nodalfieldnames,elemfieldnames,pathname,filename,part,time,binary_output)
% function [] = write_vtk_mesh(M,nodalfields,elemfields,nodalfieldnames,elemfieldnames,pathname,filename,part,time,binary_output)

if nargin<4 || isempty(nodalfieldnames)
    nodalfieldnames='solution';
end
if nargin<5 || isempty(elemfieldnames)
    elemfieldnames='solution';
end
if nargin<6 || isempty(pathname)
    pathname='.';
end
if nargin<7 || isempty(filename)
    filename='paraview';
end
if nargin<8 || isempty(part)
    part=1;
end
if nargin<9 || isempty(time)
    time=0;
end
if nargin<10 || isempty(binary_output)
    binary_output=true;
end

if ~isempty(nodalfields) && ~iscell(nodalfields)
    nodalfields={nodalfields};  
end
if ~isempty(nodalfieldnames) && ~iscell(nodalfieldnames)
    nodalfieldnames={nodalfieldnames};  
end

if ~isempty(elemfields) && ~iscell(elemfields)
    elemfields={elemfields};  
end
if ~isempty(elemfieldnames) && ~iscell(elemfieldnames)
    elemfieldnames={elemfieldnames};  
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
    types = [types repmat(type,1,getnbelem(elem))];
end

indim = getindim(M);
x = getcoord(M.node);
if indim==1
    x = [x,zeros(nnode,2)]';
elseif indim==2
    x = [x,zeros(nnode,1)]';
else
    x = x';
end
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

[filepath,filename,ext] = fileparts(filename);
if isempty(ext)
    ext = '.vtu';
end
filename = strcat(filename,'_',num2str(part),'_',num2str(time),ext);
filename = fullfile(pathname,filepath,filename);
fid = fopen(filename,'w',endian_matlab);
fprintf(fid,'<?xml version="1.0" ?>\n');

fprintf(fid,['<VTKFile type="UnstructuredGrid" version="0.1" byte_order="',endian_paraview,'">\n']);
fprintf(fid,'\t <UnstructuredGrid>\n');
fprintf(fid,'\t\t <Piece NumberOfPoints="%u" NumberOfCells="%u">\n',nnode,nelem);

if binary_output
    
    % POINT DATA
    fprintf(fid,'\t\t\t <PointData scalars="scalar"> \n');
    if ~isempty(nodalfields)
        for i=1:length(nodalfields)
            nbofcomponents = numel(nodalfields{i})/nnode;
            if nbofcomponents==2
                nodalfields{i} = reshape(nodalfields{i},nbofcomponents,nnode);
                nodalfields{i} = [nodalfields{i};zeros(1,nnode)];
                nbofcomponents = 3;
            end
            fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="%s" NumberOfComponents="%u" format="appended" offset="%u" />\n',nodalfieldnames{i},nbofcomponents,offset);
            offset=offset+4+4*numel(nodalfields{i});
        end
    end
    fprintf(fid,'\t\t\t </PointData> \n');
    
    % CELL DATA
    fprintf(fid,'\t\t\t <CellData> \n');
    if ~isempty(elemfields)
        for i=1:length(elemfields)
            nbofcomponents = numel(elemfields{i})/nelem;
            fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="%s" NumberOfComponents="%u" format="appended" offset="%u" />\n',elemfieldnames{i},nbofcomponents,offset);
            offset=offset+4+4*numel(elemfields{i});
        end
    end
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
    
    % END VTK FILE
    fprintf(fid,'\t\t </Piece>\n');
    fprintf(fid,'\t </UnstructuredGrid> \n');
    
    % APPENDED DATA
    fprintf(fid,'\t <AppendedData encoding="raw"> \n _');
    
    % NODAL FIELDS
    if ~isempty(nodalfields)
        for i=1:length(nodalfields)
            fwrite(fid,4*numel(nodalfields{i}),'uint32');
            fwrite(fid,nodalfields{i},'float32');
        end
    end
    
    % ELEM FIELDS
    if ~isempty(elemfields)
        for i=1:length(elemfields)
            fwrite(fid,4*numel(elemfields{i}),'uint32');
            fwrite(fid,elemfields{i},'float32');
        end
    end
    
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
    if ~isempty(nodalfields)
        for i=1:length(nodalfields)
            nbofcomponents = numel(nodalfields{i})/nnode;
            if nbofcomponents==2
                nodalfields{i} = reshape(nodalfields{i},nbofcomponents,nnode);
                nodalfields{i} = [nodalfields{i};zeros(1,nnode)];
                nbofcomponents = 3;
            end
            fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="%s" NumberOfComponents="%u" format="ascii">\n',nodalfieldnames{i},nbofcomponents);
            fprintf(fid,'%e \n',nodalfields{i});
            fprintf(fid,'\t\t\t\t </DataArray>\n');
        end
    end
    fprintf(fid,'\t\t\t </PointData> \n');
    
    % CELL DATA
    fprintf(fid,'\t\t\t <CellData> \n');
    if ~isempty(elemfields)
        for i=1:length(elemfields)
            nbofcomponents = numel(elemfields{i})/nelem;
            fprintf(fid,'\t\t\t\t <DataArray type="Float32" Name="%s" NumberOfComponents="%u" format="ascii">\n',elemfieldnames{i},nbofcomponents);
            fprintf(fid,'%e \n',elemfields{i});
            fprintf(fid,'\t\t\t\t </DataArray>\n');
        end
    end
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
