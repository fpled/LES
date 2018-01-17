clc
clearvars
close all

%%%%%%%%%% TecIO Setup %%%%%%%%%%%%%
if ismac
    tecplot_home = '/Applications/Tecplot 360 EX 2017 R2/Tecplot 360 EX 2017 R2.app/Contents/MacOS');
    tecio_path = strcat(tecplot_home,'/libtecio.dylib');
elseif isunix
    %tecplot_home = '/usr/local/apps/tecplot360';
    %tecio_path = strcat(tecplot_home,'/lib/libtec.so');
    tecplot_home = '/home/p/pled/Applications/Tecplot/360ex_2017r3';
    tecio_path = strcat(tecplot_home,'/bin/libtecio.so');
end
tecio_header_path = 'TECIO.h'; % uses simplified TECIO.h included in package. 
if ~libisloaded('tecio')
    [notfound,warnings]=loadlibrary(tecio_path,tecio_header_path,'alias','tecio');
end
% libfunctionsview('tecio')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varnames = {'vitesse_x1','vitesse_x2','vitesse_x3','fluide2'}; % variable names
n = length(varnames); % number of variables
p = 50; % number of time steps
N = 19; % number of samples

% for g=2.^(4:8)
for g=2.^(4:7)
    gridname = fullfile('Data',['Grid' num2str(g)]);
    m = (g+1)^3; % number of spatial points
    Y = zeros(N,m*n*p);
    for t=0:p
        filename = ['InvPhase3d_' num2str(t,'%05d00') '.szplt'];
        
        Yt = zeros(N,m*n);
        for l=1:N
            foldername = ['CasInvPhase3D-' num2str(l)];
            
            file = fullfile(gridname,foldername,filename);
            [isok,~,handle] = calllib('tecio','tecFileReaderOpen',file,[]);
            
            title = libpointer('stringPtrPtr',cell(1,1));
            [isok,~,title] = calllib('tecio','tecDataSetGetTitle',handle,title);
            
            numvars = 0;
            numzones = 0;
            [isok,~,numvars] = calllib('tecio','tecDataSetGetNumVars',handle,numvars);
            [isok,~,numzones] = calllib('tecio','tecDataSetGetNumZones',handle,numzones);
            if numvars<n
                error(['Wrong number of variables in file ' file])
            end
            if numzones>1
                error(['Wrong number of zones in file ' file])
            end
            
            Yl = zeros(m*n,1);
            for var=1:numvars
                name = libpointer('stringPtrPtr',cell(1,1));
                [isok,~,name] = calllib('tecio','tecVarGetName',handle,var,name);
                [isvar,i] = ismember(name,varnames);
                if isvar
                    zone = 1;
                    % for zone=1:numzones
                        numvals = 0;
                        [isok,~,numvals] = calllib('tecio','tecZoneVarGetNumValues',handle,zone,var,numvals);
                        if numvals~=m
                            error(['Wrong dataset for variable ' name{:} ' in file ' file])
                        end
                        
                        values = zeros(numvals,1);
                        [isok,~,values] = calllib('tecio','tecZoneVarGetFloatValues',handle,zone,var,1,numvals,values);
                        
                        Yl((0:m-1)*n+i) = values;
                    % end
                end
            end
            Yt(l,:) = Yl;
            
            calllib('tecio','tecFileReaderClose',handle) ;
        end
        Y(:,m*n*t+1:m*n*(t+1)) = Yt;
        
    end
    save(fullfile('Data',['data' num2str(g) '.mat']),'Y');
end

if libisloaded('tecio')
    unloadlibrary('tecio') 
end

