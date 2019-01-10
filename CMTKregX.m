function [] = CMTKregX(varargin)
% Affine and Warp registration of images onto a reference using CMTK's
% registrationX and WarpX.
%
% Correlation between output registration and reference image can also be
% computed (set 'Correlate' value to 1)
%
% Requirements:
%
% 1. Cygwin and CMTK should be locally installed and their directories set below.
% 2. The image channels should be named as: image1_01 for the reference
% channel, and image1_0N for the other Nth other channels.
% 3. Images are required to be .nrrd or .nii file type (gz compression is
% also supported)
% 4. Images need to have pixel size properties set.
%
% Input arguments should obey a ('Parameter', value), structure. 
% eg.,
% CMTKregX('Channels',3,'Metric',{'nmi','ncc'},'Auto_Multi_Levels',[0,1])
%
% Parameters not addressed will be set to default values (see below).
%
% v2 - adds iteration through Init parameter
% v3 - commands can now deal with strings that have spaces in it. Added
% ability to use .nii extensions too.
% v4 - warpx elastic transformation now supported.
% 170913 - ignore registration if already exists and addition of float and
% fixed images as possible input
% 171122 - Windows 10 bash now supported

w10bash = 0;    % Using Windows 10 bash?
if w10bash
    bashdir = 'bash';
    cmtkdir = '/usr/lib/cmtk/bin';  % CMTK linux instalation directory
    mntcmd = 'sudo mount -t drvfs ''\\\128.40.155.187\data2'' /mnt/data2';  % server path to mount in bash
else
    bashdir = 'C:\cygwin64\bin\bash';   % Cygwin \bin\bash directory
    cmtkdir = 'C:/Fiji.app/bin/cmtk';   % CMTK instalation directory    
end

%%  Set-up registration parameters

parameters.Registration = 1;    % Compute linear registration? (1/0)
parameters.Warp = 0;    % Compute warp registration? (1/0)
parameters.Reformat = 1;    % Appy transformation to images? (1/0)
parameters.Channels = 2;    % Number of channels
parameters.FloatDir = [];   % Floating images directory
parameters.FixedFile = [];  % Reference filepath

% Default parameter values for registrationx
parameters.Auto_Multi_Levels = 0;
parameters.Max_Stepsize = 250;
parameters.Min_Stepsize = 0.1;
parameters.Stepfactor = 0.1; % 0-1
parameters.Metric = {'nmi'}; % ncc/nmi/mi/cr/rms/msd
parameters.Coarsest = -1;
parameters.Dofs = '--dofs 6 --dofs 9'; % 12 = shearing
parameters.Init = {'com'}; % none/fov/com/pax

% Default parameter values for warpx
parameters.Grid = 80;
parameters.Energy = 1e-1;
parameters.Refine = 3;
parameters.Fastness = 'fast';
parameters.WMetric = 'ncc';
parameters.WCoarsest = 8;
parameters.WExploration = 10;
parameters.WAccuracy = -1;

if nargin > 0
    parnames = fieldnames(parameters);
    
    % look for right argument structure
    if mod(nargin,2) ~= 0
        error('Arguments must occur in name-value pairs.');
    end
    
    for k = 1:2:nargin
        arg = varargin{k};  % Parameter name
        if ~ischar(arg)
            error('Expected argument %d to be a string parameter name.', k);
        end
        if any(strcmp(arg,parnames))
            parameters.(arg) = varargin{k+1};
        else
            error('Unrecognized parameter name ''%s''.', arg);
        end
    end
else
    disp('Using default registration parameters...')
end

% Show parameters
parameters

%% Start the registration

% Image folders to use
if isempty(parameters.FloatDir)
    folders = uipickfiles('FilterSpec','\\128.40.168.141\bdata2\Dammy\PAGFP', ...
        'Prompt','Select folders for registration');
else
    folders = {parameters.FloatDir};
end

% Reference image to use
if isempty(parameters.FixedFile)
    [refname, refdir] = uigetfile('*.n*','Select Reference brain file','C:\mapping\References');
else
    refstr = strsplit(parameters.FixedFile,filesep);
    if isempty(refstr{1})
        refdir = [filesep filesep fullfile(refstr{1:end-1})];
    else
        refdir = fullfile(refstr{1:end-1});
    end
    refname = refstr{end};
end

reffile = strrep(fullfile(refdir,refname),filesep,'/');
refsplt = strsplit(reffile,'.');
if length(refsplt) > 2
    refext = refsplt{2}; else; refext = refsplt{end}; end

for f = 1:length(folders)
    
    datadir = folders{f};
    regdir = strrep(datadir,filesep,'/');        % for cygwin
    
    fltfiles1 = dir([datadir '\' '*_01.*']);
    
    if size(fltfiles1,1) > 1
        sprintf('More that one "01" channel present... Cannot proceed for: \n%s', datadir)
    else
        fltnames = strsplit(fltfiles1.name,'.');
        fltname = fltnames{1};
        if length(fltnames) == 3
            fltext = [fltnames{2} '.' fltnames{3}];
        else
            fltext = fltnames{end};
        end
        fltfile = [regdir '/' fltname];
        
        k = 1;
        for aml = 1:length(parameters.Auto_Multi_Levels)
            for mxs = 1:length(parameters.Max_Stepsize)
                for mns = 1:length(parameters.Min_Stepsize)
                    for sf = 1:length(parameters.Stepfactor)
                        for mtr = 1:length(parameters.Metric)
                            for cor = 1:length(parameters.Coarsest)
                                for i = 1:length(parameters.Init)
                                    for grd = 1:length(parameters.Grid)
                                        for ene = 1:length(parameters.Energy)
                                            for refn = 1:length(parameters.Refine)
                                                for wacc = 1:length(parameters.WAccuracy)
                                                    for wexp = 1:length(parameters.WExploration)
                                                        for wcor = 1:length(parameters.WCoarsest)
                                                            
                                                            % Unique file name
                                                            tname_affine = sprintf('aml%g_mxs%g_mns%g_sf%g_%s_c%g_i%s', ...
                                                                parameters.Auto_Multi_Levels(aml), parameters.Max_Stepsize(mxs), ...
                                                                parameters.Min_Stepsize(mns), parameters.Stepfactor(sf), ...
                                                                parameters.Metric{mtr}, parameters.Coarsest(cor), parameters.Init{i});
                                                            
                                                            tname_warp = sprintf('grd%g_ene%g_ref%g_coar%g_exp%g_acc%g', ...
                                                                parameters.Grid(grd), parameters.Energy(ene), parameters.Refine(refn), ...
                                                                parameters.WCoarsest(wcor), parameters.WExploration(wexp), parameters.WAccuracy(wacc));
                                                            
                                                            if parameters.Registration == 1 && parameters.Warp == 1
                                                                tname = [tname_affine '_' tname_warp];
                                                            elseif parameters.Registration == 1 && parameters.Warp == 0
                                                                tname = tname_affine;
                                                            elseif parameters.Registration == 0 && parameters.Warp == 1
                                                                tname = tname_warp;
                                                            end
                                                            
                                                            outdir = fullfile(datadir,tname);
                                                            if ~exist(outdir,'dir')
                                                                mkdir(outdir);
                                                            else
                                                                % Skip if file dir already exists
                                                                if length(dir(fullfile(datadir,tname))) > 2
                                                                    continue
                                                                end
                                                            end
                                                            
                                                            routfile = sprintf('./%s/%s_affine.xform',tname,fltname);
                                                            matfile = sprintf('./%s/%s_affine_mat.csv',tname,fltname);
                                                            woutfile = sprintf('./%s/%s_warp.xform',tname,fltname);
                                                            
                                                            % Linear registration
                                                            if parameters.Registration == 1
                                                                
                                                                fprintf('Registering %s (%d)\n',datadir,(k))
                                                                
                                                                rregcmd = sprintf(['-v %s --auto-multi-levels %d --max-stepsize %d --min-stepsize %d --symmetric ' ...
                                                                    '--stepfactor %d --%s --coarsest %d --%s -o ''%s'' --write-matrix ''%s'' ''%s'' ''%s.%s'''], ...
                                                                    parameters.Dofs, ...
                                                                    parameters.Auto_Multi_Levels(aml), ...
                                                                    parameters.Max_Stepsize(mxs), ...
                                                                    parameters.Min_Stepsize(mns), ...
                                                                    parameters.Stepfactor(sf), ...
                                                                    parameters.Metric{mtr}, ...
                                                                    parameters.Coarsest(cor), ...
                                                                    parameters.Init{i}, ...
                                                                    routfile, matfile, reffile, fltfile, fltext);
                                                                
                                                                if w10bash
                                                                    syscmd = sprintf('%s --login -c "%s; cd ''%s''; %s/registrationx %s &"', ...
                                                                        bashdir,mntcmd,regdir,cmtkdir,rregcmd);
                                                                    syscmd = strrep(syscmd,filesep,'/');    % correct for unix slash
                                                                    syscmd = strrep(syscmd,'///128.40.155.187/data2','\\\128.40.155.187\data2');
                                                                    
                                                                    % Correct for bash mounting
                                                                    driveltrpos = strfind(syscmd,':/');
                                                                    if ~isempty(driveltrpos)
                                                                        driveltr = syscmd(driveltrpos(1)-1);    % mount drive letter
                                                                        syscmd = strrep(syscmd,sprintf('%s:/',driveltr), ...
                                                                            sprintf('/mnt/%s/',lower(driveltr)));
                                                                    end
                                                                    
                                                                    % Correct for path in server
                                                                    syscmd = strrep(syscmd,'//128.40.155.187/data2/','/mnt/data2/');
                                                                else
                                                                    syscmd = sprintf('%s --login -c "cd ''%s''; %s/registrationx %s"', ...
                                                                        bashdir,regdir,cmtkdir,rregcmd);
                                                                end
                                                                
                                                                % Run command
                                                                system(syscmd,'-echo');
                                                                
                                                                % Convert to ITK format 
                                                                itkoutfile = strrep(routfile,'.xform','_itk.tfm');
                                                                syscmd = sprintf(['%s --login -c "cd ''%s''; ' ...
                                                                    '%s/xform2itk --fixed-space ''LPS'' --moving-space ''LPS'' %s %s"'], ...
                                                                    bashdir,regdir,cmtkdir,routfile,itkoutfile);
                                                                system(syscmd,'-echo');
                                                                
                                                                if w10bash
                                                                    % Check for last loop finish
                                                                    finished = 0;
                                                                    while ~finished
                                                                        file = dir(fullfile(outdir,sprintf('%s_affine_mat.csv',fltname)));
                                                                        if ~isempty(file)
                                                                            finished = 1;
                                                                        else
                                                                            pause(10); % check every 10s
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                            
                                                            % Warping
                                                            if parameters.Warp == 1
                                                                
                                                                fprintf('Warping %s (%d)\n',datadir,(k))
                                                                
                                                                wregcmd = sprintf(['-v --%s --grid-spacing %d -e %d --energy-weight %d ' ...
                                                                    '--%s --coarsest %d --accuracy %d --refine %d -o ''%s'' ''%s'''], ...
                                                                    parameters.Fastness, ...
                                                                    parameters.Grid(grd), ...
                                                                    parameters.WExploration(wexp), ...
                                                                    parameters.Energy(ene), ...
                                                                    parameters.WMetric, ...
                                                                    parameters.WCoarsest(wcor), ...
                                                                    parameters.WAccuracy(wacc), ...
                                                                    parameters.Refine(refn), ...
                                                                    woutfile,routfile);
                                                                
                                                                if w10bash
                                                                    syscmd = sprintf('%s --login -c "%s; cd ''%s''; %s/warp %s" &', ...
                                                                        bashdir,mntcmd,regdir,cmtkdir,wregcmd);
                                                                    syscmd = strrep(syscmd,filesep,'/');    % correct for unix slash
                                                                    syscmd = strrep(syscmd,'///128.40.155.187/data2','\\\128.40.155.187\data2');
                                                                    
                                                                    % Correct for bash mounting
                                                                    driveltrpos = strfind(syscmd,':/');
                                                                    if ~isempty(driveltrpos)
                                                                        driveltr = syscmd(driveltrpos(1)-1);    % mount drive letter
                                                                        syscmd = strrep(syscmd,sprintf('%s:/',driveltr), ...
                                                                            sprintf('/mnt/%s/',lower(driveltr)));
                                                                    end
                                                                    % Correct for path in server
                                                                    syscmd = strrep(syscmd,'//128.40.155.187/data2/','/mnt/data2/');
                                                                else
                                                                    syscmd = sprintf('%s --login -c "cd ''%s''; %s/warp %s"', ...
                                                                        bashdir,regdir,cmtkdir,wregcmd);
                                                                end
                                                                
                                                                % Run command
                                                                system(syscmd,'-echo');
                                                                
                                                                if w10bash
                                                                    Check for last loop finish
                                                                    finished = 0;
                                                                    while ~finished
                                                                        file = dir(fullfile(outdir,sprintf('%s_affine_mat.csv',fltname)));
                                                                        if ~isempty(file)
                                                                            finished = 1;
                                                                        else
                                                                            pause(10); % check every 10s
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                            
                                                            % Loop through all channels
                                                            for nc = 1:parameters.Channels
                                                                
                                                                routreff = sprintf('./%s/%s',tname, ...
                                                                    strrep(fltname, '_01', sprintf('_0%d_affine.%s',nc,fltext)));
                                                                woutreff = sprintf('./%s/%s',tname, ...
                                                                    strrep(fltname, '_01', sprintf('_0%d_warp.%s',nc,fltext)));
                                                                inreff = strrep(fltname, '_01', sprintf('_0%d.%s',nc,fltext));
                                                                
                                                                % Reformat
                                                                if parameters.Reformat == 1
                                                                    
                                                                    % Reformat
                                                                    % affine
                                                                    sprintf('Reformating affine for Channel %d ...',nc)
                                                                    rrregcmd = sprintf('-v -o ''%s'' --floating ''%s'' ''%s'' ''%s''', ...
                                                                        routreff, inreff, reffile, routfile);
                                                                    
                                                                    if w10bash
                                                                        syscmd = sprintf('%s --login -c "%s; cd ''%s''; %s/reformatx %s" &', ...
                                                                            bashdir,mntcmd,regdir,cmtkdir,rrregcmd);
                                                                        syscmd = strrep(syscmd,filesep,'/');    % correct for unix slash
                                                                        syscmd = strrep(syscmd,'///128.40.155.187/data2','\\\128.40.155.187\data2');
                                                                        
                                                                        % Correct for bash mounting
                                                                        driveltrpos = strfind(syscmd,':/');
                                                                        if ~isempty(driveltrpos)
                                                                            driveltr = syscmd(driveltrpos(1)-1);    % mount drive letter
                                                                            syscmd = strrep(syscmd,sprintf('%s:/',driveltr), ...
                                                                                sprintf('/mnt/%s/',lower(driveltr)));
                                                                        end
                                                                        % Correct for path in server
                                                                        syscmd = strrep(syscmd,'//128.40.155.187/data2/','/mnt/data2/');
                                                                    else
                                                                        syscmd = sprintf('%s --login -c "cd ''%s''; %s/reformatx %s"', ...
                                                                            bashdir,regdir,cmtkdir,rrregcmd);
                                                                    end
                                                                    
                                                                    system(syscmd,'-echo');
                                                                    
                                                                    if w10bash
                                                                        % Check for last loop finish
                                                                        finished = 0;
                                                                        while ~finished
                                                                            file = dir(fullfile(outdir,sprintf('%s_affine_mat.csv',fltname)));
                                                                            if ~isempty(file)
                                                                                finished = 1;
                                                                            else
                                                                                pause(10); % check every 10s
                                                                            end
                                                                        end
                                                                    end
                                                                    
                                                                    % Reformat
                                                                    % warp
                                                                    if parameters.Warp == 1
                                                                        sprintf('Reformating warp for Channel %d ...',nc)
                                                                        wrregcmd = sprintf('-v -o ''%s'' --floating ''%s'' ''%s'' ''%s''', ...
                                                                            woutreff, inreff, reffile, woutfile);
                                                                        
                                                                        if w10bash
                                                                            syscmd = sprintf('%s --login -c "%s; cd ''%s''; %s/reformatx %s" &', ...
                                                                                bashdir,mntcmd,regdir,cmtkdir,wrregcmd);
                                                                            syscmd = strrep(syscmd,filesep,'/');    % correct for unix slash
                                                                            syscmd = strrep(syscmd,'///128.40.155.187/data2','\\\128.40.155.187\data2');
                                                                            
                                                                            % Correct for bash mounting
                                                                            driveltrpos = strfind(syscmd,':/');
                                                                            if ~isempty(driveltrpos)
                                                                                driveltr = syscmd(driveltrpos(1)-1);    % mount drive letter
                                                                                syscmd = strrep(syscmd,sprintf('%s:/',driveltr), ...
                                                                                    sprintf('/mnt/%s/',lower(driveltr)));
                                                                            end
                                                                            % Correct for path in server
                                                                            syscmd = strrep(syscmd,'//128.40.155.187/data2/','/mnt/data2/');
                                                                        else
                                                                            syscmd = sprintf('%s --login -c "cd ''%s''; %s/reformatx %s"', ...
                                                                                bashdir,regdir,cmtkdir,wrregcmd);
                                                                        end
                                                                        
                                                                        system(syscmd,'-echo');
                                                                        
                                                                        if w10bash
                                                                            % Check for last loop finish
                                                                            finished = 0;
                                                                            while ~finished
                                                                                file = dir(fullfile(outdir,sprintf('%s_affine_mat.csv',fltname)));
                                                                                if ~isempty(file)
                                                                                    finished = 1;
                                                                                else
                                                                                    pause(10); % check every 10s
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                            k = k+1;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end