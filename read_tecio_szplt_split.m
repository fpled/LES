clc
clearvars
close all

%%%%%%%%%% TecIO Setup %%%%%%%%%%%%%
if ismac
    tecplot_home = '/Applications/Tecplot 360 EX 2017 R2/Tecplot 360 EX 2017 R2.app/Contents/MacOS';
    tecio_path = strcat(tecplot_home,'/libtecio.dylib');
elseif isunix
    % tecplot_home = '/usr/local/apps/tecplot360';
    % tecio_path = strcat(tecplot_home,'/lib/libtec.so');
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
N = 40; % number of samples

pathname = fileparts(mfilename('fullpath'));
% pathname = '/mnt/tcm13/SV_FP/';

for g=2.^(7:8)
    gridname = ['Grid' num2str(g)];
    disp(gridname)
    
    m = (g+1)^3; % number of spatial points
    for t=0:p
        time = ['Time ' num2str(t)];
        disp(time)
        
        filename = ['InvPhase3d_' num2str(t,'%05d00') '.szplt'];
        
        % Yt = zeros(N,n*m);
        Yt = zeros(N,n,m);
        for l=1:N
            foldername = ['Cas-' num2str(l)];
            disp(foldername)
            
            file = fullfile(pathname,gridname,foldername,'validation',filename);
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
            
            % Yl = zeros(m*n,1);
            Yl = zeros(n,m);
            for var=1:numvars
                name = libpointer('stringPtrPtr',cell(1,1));
                [isok,~,name] = calllib('tecio','tecVarGetName',handle,var,name);
                [isvar,i] = ismember(name,varnames);
                if isvar
                    zone = 1;
                    % for zone=1:numzones
%                         type = 0;
%                         [isok,~,type] = calllib('tecio','tecZoneGetType',handle,zone,type);
%                         
%                         % 0 - Ordered, 1- FE line, ... 5- FE Brick...
%                         if type == 0
%                             I = 0; J = 0; K = 0;
%                             [isok,~,I,J,K] = calllib('tecio','tecZoneGetIJK',handle,zone,I,J,K);
%                         end
                        
                        numvals = 0;
                        [isok,~,numvals] = calllib('tecio','tecZoneVarGetNumValues',handle,zone,var,numvals);
                        if numvals~=m
                            error(['Wrong dataset for variable ' name{:} ' in file ' file])
                        end
                        
                        values = zeros(numvals,1);
                        [isok,~,values] = calllib('tecio','tecZoneVarGetFloatValues',handle,zone,var,1,numvals,values);
                        
                        % Yl((0:m-1)*n+i) = values;
                        Yl(i,:) = values;
                    % end
                end
            end
            % Yt(l,:) = Yl;
            Yt(l,:,:) = Yl;
            
            calllib('tecio','tecFileReaderClose',handle);
        end
        save(fullfile(pathname,gridname,['data_t' num2str(t) '.mat']),'Yt');
        
    end
    fprintf('\n');
    save(fullfile(pathname,gridname,'data.mat'),'N','n','m','p');
end

if libisloaded('tecio')
    unloadlibrary('tecio')
end

