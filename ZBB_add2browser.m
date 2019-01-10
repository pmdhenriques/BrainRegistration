function [] = ZBB_add2browser()
% Adds ZBB registered 3D image stacks to the ZBB browser
%
% Pedro Henriques, Sep 2017

zbbimsize = [1030,616,420]; % ZBB image size
zbbimdir = '\\128.40.168.141\bdata2\Registration\ZBB_browser_v2\stack_images'; % ZBB stack_images directory

dirchoice = questdlg('Add new directory?', ...
    'New directory', ...
    'List','New','List');
switch dirchoice
    case 'List'
        ZBBfolders = dir(zbbimdir);
        ZBBfolders = {ZBBfolders([ZBBfolders.isdir] == 1).name};
        ZBBfolders = ZBBfolders(3:end);
        ok = 0;
        while ~ok
            [diridx,ok] = listdlg('PromptString','Select a mask:',...
                'SelectionMode','single',...
                'ListString',ZBBfolders(:));
        end
        outdir = fullfile(zbbimdir,ZBBfolders{diridx});
    case 'New'
        dirstr = inputdlg('Plese write the name');
        dirstr = dirstr{1};
        mkdir(fullfile(zbbimdir,dirstr));
        outdir = fullfile(zbbimdir,dirstr);
end

%%

[filename, filepath] = uigetfile(fullfile( ...
    '\\128.40.155.187\data2\Bianco_lab\Pedro\NI project\NI_Electroporations', ...
    '*.nii*'), ...
    'Select nii file');

fileprefix = strsplit(filename,'.nii');
fileprefix = fileprefix{1};

h = waitbar(0,'Loading');
im = load_untouch_nii(fullfile(filepath,filename));
im = im.img;

if isequal(size(im),zbbimsize)
    
    im(im<0) = 0;
    im = rot90(uint8(scaledata(im,0,255)));
    
    %% Write TIF
    
    tf = Tiff(fullfile(outdir,[fileprefix '.tif']), 'w');
    tags.Photometric = Tiff.Photometric.MinIsBlack;
    tags.Compression = 5;
    tags.ImageLength = size(im, 1);
    tags.ImageWidth = size(im, 2);
    tags.SampleFormat = Tiff.SampleFormat.UInt;
    tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tags.Orientation = Tiff.Orientation.TopLeft;
    tags.BitsPerSample = 8;
    
    nfr = size(im,3);
    for fr = 1:nfr
        waitbar(fr/nfr,h,'Saving TIF');
        tf.setTag(tags)
        tf.write(im(:,:,fr));
        if fr ~= nfr
            tf.writeDirectory()
        end
    end
    tf.close()
    close(h);
else
    errordlg('Image dimentions don''t match');
end
end