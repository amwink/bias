function fegui;
%
% graphical user interface to the fastECM program
%
%
%
% (c) Alle Meije Wink -- 16/03/2012
%     a.m.winkATgmail.com
%
% If you use this program in your research, remember to cite this paper
% in the journal "Brain Connectivity":
%
% Alle Meije Wink, Jan C de Munck, Ysbrand D van der Werf, Odile A van den heuvel, Frederik Barkhof
%         "Fast eigenvector centrality mapping of voxel-wise connectivity in functional MRI:
%          implementation, validation and interpretation."
%         Brain Connectivity 2012, Vol. 2 No. 5, pages 265-274
% URL:    http://online.liebertpub.com/doi/abs/10.1089/brain.2012.0087
% or:     http://dare.ubvu.vu.nl/handle/1871/48750
%
%
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% build the window & gui objects
%

% delete existing guis and make current one

fegui = findobj('tag', 'fegui');
delete(fegui);

fegui = figure('tag', 'fegui', ...
    'menubar', 'none', ...
    'position', [100 100 750 450], ...
    'name', 'fast ECM', 'numbertitle', 'off');

% make 'Files' and 'Settings' tabs

warning off
hTabGroup = uitabgroup;
ftab      = uitab(hTabGroup, 'title', '                Files                ');
stab      = uitab(hTabGroup, 'title', '               Settings              ');
warning on

h.exts = {  '*.nii;*.nii.gz;*.img;*.img.gz' 'all image types'              ...
    '*.nii'                         'nifti files'                  ...
    '*.nii.gz'                      'compressed nifti files'       ...
    '*.img'                         'analyze files'                ...
    '*.img.gz'                      'compressed analyze files'     ...
    '*.txt'                         'text files with file names'   };
h.exts = reshape(h.exts, 2, 6)';

% add components

h.gb = uicontrol (fegui, ...
    'style',                'pushbutton', ...
    'string',               'estimate fECM', ...
    'tag',                  'gobutton', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tooltipstring',        'start computing the fECM of each input in the current selection', ...
    'Position',             [.05 .85 .4 .05 ]);

h.sb = uicontrol (ftab, ...
    'style',                'pushbutton', ...
    'string',               'add files', ...
    'tag',                  'selectbutton', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tooltipstring',        'add files to the selection of fECM inputs', ...
    'Position',             [.05 .8 .4 .05 ]);

h.fc = uicontrol (ftab, ...
    'style',                'text', ...
    'string',               'selected files:', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'thefilescaption', ...
    'horizontalalignment',  'left', ...
    'tooltipstring',        'current selection of files whose fECM will be computed', ...
    'position',             [.55 .9 .4 .05]);

h.fv = uicontrol (ftab, ...
    'style',                'listbox', ...
    'string',               '', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'thefilesvalue', ...
    'horizontalalignment',  'left', ...
    'tooltipstring',        sprintf('file list\nleft click: save to file\nright click: clear'), ...
    'position',             [.55 .10 .4 .75]);

h.mb = uicontrol (ftab, ...
    'style',                'pushbutton', ...
    'string',               'use mask', ...
    'tag',                  'maskbutton', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tooltipstring',        sprintf('use specified brain mask\n(default: none - use nonzero voxels)'), ...
    'Position',             [.05 .2 .4 .05 ]);

h.mv = uicontrol (ftab, ...
    'style',                'text', ...
    'string',               ' no mask', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'maskvalue', ...
    'backgroundcolor',      [.8 .8 .8], ...
    'horizontalalignment',  'left', ...
    'tooltipstring',        sprintf('mask file\nright click: clear'), ...
    'position',             [.05 .11 .4 .04]);

h.ab = uicontrol (ftab, ...
    'style',                'pushbutton', ...
    'string',               'use atlas', ...
    'tag',                  'atlasbutton', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tooltipstring',        sprintf('use atlas for regional analysis\n(default: none - use nonzero voxels)'), ...
    'Position',             [.05 .4 .4 .05 ]);

h.av = uicontrol (ftab, ...
    'style',                'text', ...
    'string',               ' no atlas', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'atlasvalue', ...
    'backgroundcolor',      [.8 .8 .8], ...
    'horizontalalignment',  'left', ...
    'tooltipstring',        sprintf('atlas file\nright click: clear'), ...
    'position',             [.05 .31 .4 .04]);

h.pv = uicontrol (stab, ...
    'style',                'checkbox', ...
    'value',                0, ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'powervalue', ...
    'horizontalalignment',  'left', ...
    'position',             [.05 .61 .05 .05]);

h.pc = uicontrol (stab, ...
    'style',                'text', ...
    'string',               'write node power map *degCM.nii*', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'powercaption', ...
    'horizontalalignment',  'left', ...
    'tooltipstring',        sprintf('produce maps of node power:\nsum of connection strengths per node'), ...
    'position',             [.1 .6 .4 .05]);

h.rv = uicontrol (stab, ...
    'style',                'checkbox', ...
    'value',                0, ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'rankvalue', ...
    'horizontalalignment',  'left', ...
    'position',             [.05 .51 .05 .05]);

h.rc = uicontrol (stab, ...
    'style',                'text', ...
    'string',               'write uniformly distributed map *rankECM.nii*', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'rankcaption', ...
    'horizontalalignment',  'left', ...
    'tooltipstring',        sprintf('produce fECM rank maps, with\na uniform distribution on ]0, 1['), ...
    'position',             [.1 .5 .4 .05]);

h.nv = uicontrol (stab, ...
    'style',                'checkbox', ...
    'value',                0, ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'normvalue', ...
    'horizontalalignment',  'left', ...
    'position',             [.05 .41 .05 .05]);

h.nc = uicontrol (stab, ...
    'style',                'text', ...
    'string',               'write normally distributed map *normECM.nii*', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'normcaption', ...
    'horizontalalignment',  'left', ...
    'tooltipstring',        sprintf('compute fECM maps with N(0, 1)\ndistributed values, via fECM ranks'), ...
    'position',             [.1 .4 .4 .05]);

h.ic = uicontrol (stab, ...
    'style',                'text', ...
    'string',               'max. iterations', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'itercaption', ...
    'horizontalalignment',  'left', ...
    'tooltipstring',        sprintf('maximum number of iterations\n(in the case of no convergence)\nfor the fECM algorithm'), ...
    'position',             [.05 .305 .45 .05]);

h.iv = uicontrol (stab, ...
    'style',                'slider', ...
    'min',                  0, ...
    'max',                  100, ...
    'value',                50, ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'itervalue', ...
    'horizontalalignment',  'left', ...
    'position',             [.2 .31 .2 .05]);

h.ie = uicontrol (stab, ...
    'style',                'edit', ...
    'visible',              'on', ...
    'string',               '50', ...
    'units',                'normalized', ...
    'tag',                  'iteredit', ...
    'horizontalalignment',  'left', ...
    'position',             [.425 .31 .05 .05]);

h.dc = uicontrol (stab, ...
    'style',                'text', ...
    'string',               '# dynamics', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'dyncaption', ...
    'horizontalalignment',  'left', ...
    'tooltipstring',        sprintf('number of dynamic ECMs (1 dynamic uses the whole time series, otherwise sliding window)'), ...
    'position',             [.3 .2 .45 .05]);

h.de = uicontrol (stab, ...
    'style',                'edit', ...
    'visible',              'on', ...
    'string',               '1', ...
    'units',                'normalized', ...
    'tag',                  'dynedit', ...
    'horizontalalignment',  'left', ...
    'position',             [.425 .21 .05 .05]);

h.wc = uicontrol (stab, ...
    'style',                'text', ...
    'string',               'use whole matrix', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'powercaption', ...
    'horizontalalignment',  'left', ...
    'tooltipstring',        sprintf('compute connections matrix instead of using fast recipe\nWARNING this is a CPU and memory intensive option!'), ...
    'position',             [.25 .1 .2 .05]);

h.wv = uicontrol (stab, ...
    'style',                'checkbox', ...
    'value',                0, ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'powervalue', ...
    'horizontalalignment',  'left', ...
    'position',             [.425 .11 .05 .05]);

h.oc = uicontrol (stab, ...
    'style',                'text', ...
    'string',               'output log:', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'thefilescaption', ...
    'horizontalalignment',  'left', ...
    'position',             [.55 .9 .4 .05]);

h.ov = uicontrol (stab, ...
    'style',                'listbox', ...
    'string',               '<empty>', ...
    'visible',              'on', ...
    'units',                'normalized', ...
    'tag',                  'theoutputvalue', ...
    'horizontalalignment',  'left', ...
    'tooltipstring',        sprintf('output log\n left click: save to file\nright click: clear'), ...
    'position',             [.55 .10 .4 .75]);

% mask and atlas dirs and working dir

h.md = uicontrol(ftab, 'style', 'text', 'string', '',  'visible', 'off');
h.ad = uicontrol(ftab, 'style', 'text', 'string', '',  'visible', 'off');
h.dv = uicontrol(ftab, 'style', 'text', 'string', pwd, 'visible', 'off');

set(h.gb, 'callback', {@doecm, h});
set(h.sb, 'callback', {@getfiles, h});
set(h.mb, 'callback', {@getmask, h});
set(h.ab, 'callback', {@getatlas, h});
set(h.fv, 'callback', {@savfiles, h});
set(h.ie, 'callback', {@upd_slider, h});
set(h.iv, 'callback', {@upd_edit, h});
set(h.ov, 'callback', {@sav, h});

set(h.mv, 'buttondownfcn', @(s, e)clrmask(h), 'value', 0);
set(h.av, 'buttondownfcn', @(s, e)clratlas(h), 'value', 0);
set(h.fv, 'buttondownfcn', @(s, e)clrfiles(h));
set(h.ov, 'buttondownfcn', @(s, e)clr(h));

% check if nifti support exists otherwise add fegui-supplied version
if (exist('load_untouch_nii')~=2)
    
    fprintf('Software for reading NifTI images not found.        \n');
    fprintf('Using Tools for Nifti/Analyze for NifTI file I / O. \n');
    fprintf('  www.mathworks.com/matlabcentral/fileexchange/8797 \n');
    npath=[fileparts(which('fegui.m')) filesep 'tools4nifti'];
    addpath(npath);                   

end % if exist

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% execute fastECM on the selected files with selected options
%

function doecm(callingObject, event, h);

fnames = get(h.fv, 'string');

if ~strcmp(fnames, '')
    
    opts.maxiter  = get(h.iv, 'value');
    opts.dynamics = str2num(get(h.de, 'string'));
    opts.rankmap  = get(h.rv, 'value');
    opts.normmap  = get(h.nv, 'value');
    opts.degmap   = get(h.pv, 'value');
    opts.wholemat = get(h.wv, 'value');

    if(get(h.mv, 'value'))      
      opts.maskfile  = [get(h.md, 'string') strtrim(get(h.mv, 'string'))];      
    end % if
    
    if(get(h.av, 'value'))
      opts.atlasfile = [get(h.ad, 'string') strtrim(get(h.av, 'string'))];
    end % if
    
    set(h.ov, 'string', '');
    set(h.ov, 'value', 1);
         
    drawnow;
    
    for i = 1:length(fnames)
        
        set(h.gb, 'string', 'see log on ''Settings'' tab', 'enable', 'off')
 	
	opts.inputfile=fnames{i};
	
	cmd=textscan(evalc('opts'),'%s');	
        set(h.ov, 'string', [get(h.ov, 'string');{['running with [ ' ...
		   sprintf('%s ',cmd{:}{:}) ' ]']}]);
        set(h.ov, 'listboxtop', length(get(h.ov, 'string')))
        
        drawnow;
        cmdlog = evalc('fastECM(opts);');
        cmdlog = strread(cmdlog, '%s', 'delimiter', '\n');
        set(h.ov, 'string', [get(h.ov, 'string');cmdlog]);
        
        drawnow;
        
    end % for i
    
    set(h.gb, 'string', 'estimate fECM', 'enable', 'on')
    
end % if strcmp

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% select more files to be analysed
%

function getfiles(callingObject, event, h);

[fnames, dname] = uigetfile(  h.exts, ...
    'select file(s) to open', ...
    'multiselect', 'on', ...
    get(h.dv, 'string') );

if (dname)
    
    if ~iscell(fnames)
        fnames = {fnames};
    end % if
    
    istxt = (strcmp(fnames{1}((end-3):end), '.txt'));
    
    if (istxt)
        
        txtinputs = [fnames];
        fnames = [];
        
        for i = 1:length(txtinputs)
            
            tst = 0;
            lines = textread([dname filesep txtinputs{i}], '%s');
            
            for l = 1:length(lines)
                
                line = lines{l};
                
                try
                    tst = load_untouch_header_only(line);
                catch,
                    fprintf('file %s could not be opened\n', line);
                end % try
                
                if(isstruct(tst))
                    
                    fnamelist = get(h.fv, 'string');
                    fnamelist{end+1} = line;
                    set(h.fv, 'string', fnamelist, 'value', 1);
                    
                    drawnow;
                    
                end % if
                
            end % for l
            
        end % for i
        
    else
        
        for i = 1:length(fnames)
            
            tst = 0;
            line = [dname fnames{i}];
            
            try
                tst = load_untouch_header_only(line);
            catch
                fprintf('file %s could not be opened\n', line);
            end
            
            if(isstruct(tst))
                
                fnamelist = get(h.fv, 'string');
                fnamelist{end+1} = line;
                set(h.fv, 'string', fnamelist, 'value', 1);
                
                drawnow;
                
            end
            
        end % for i
        
    end % if istxt
    
    set(h.dv, 'string', dname);
    
end %if dname

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% save the selection of files
%

function savfiles(callingObject, event, h);

[fname, dname] = uiputfile(   {'*.txt', 'txt files'}, ...
    'file to save to', ...
    [get(h.dv, 'string') filesep 'fegui.txt']);

if (fname)
    
    fid = fopen([dname filesep fname], 'w');
    s = get(h.fv, 'string');
    
    if ~iscell(s)
        s = {s};
    end % if
    
    for i = 1:length(s)
        fprintf(fid, '%s\n', s{i});
    end % for
    
    fclose (fid);
    
end % if fname

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% select mask
%

function getmask(callingObject, event, h);

[fnames, dname] = uigetfile(  h.exts, ...
    'select file(s) to open', ...
    'multiselect', 'off', ...
    get(h.dv, 'string') );
if (dname)
    
    line = [dname fnames];
    
    try
        tst = load_untouch_header_only(line);
    catch
        fprintf('file %s could not be opened\n', line);
    end % try
    
    if(isstruct(tst))
        
        set(h.mv, 'string', [dname filesep fnames]);
        [a b c] = fileparts(get(h.mv, 'string'));
        set(h.mv, 'string', [' ' b c], 'value', 1);
        set(h.md, 'string', a);
        
        drawnow;
        
    end % if
    
end % if dname

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% select atlas
%

function getatlas(callingObject, event, h);

[fnames, dname] = uigetfile(h.exts, ...
    'select file(s) to open', ...
    'multiselect', 'off', ...
    get(h.dv, 'string') );
if (dname)
    
    line = [dname fnames];
    
    try
        tst = load_untouch_header_only(line);
    catch
        fprintf('file %s could not be opened\n', line);
    end % try
    
    if(isstruct(tst))
        
        set(h.av, 'string', [dname filesep fnames]);
        [a b c] = fileparts(get(h.av, 'string'));
        set(h.av, 'string', [' ' b c], 'value', 1 );
        set(h.ad, 'string', a);
        
        drawnow;
        
    end % if
    
end % if dname

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% change max iterations slider value after changing the edit box
%

function upd_slider(callingObject, event, h);

if( prod(size(str2num(get(h.ie, 'string')))) == 1 )
    
    iv = max(1, min(100, str2num(get(h.ie, 'string'))));
    set(h.ie, 'string', num2str(iv));
    set(h.iv, 'value', iv);
    
else
    
    upd_edit(h.iv, 0, h)
    
end % if

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% change max iterations edit box value after changing the slider
%

function upd_edit(callingObject, event, h);

set(h.iv, 'value', round(get(h.iv, 'value')));
ie = get(h.iv, 'value');
set(h.ie, 'string', num2str(ie));

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% save the output log
%

function sav(callingObject, event, h);

[fname, dname] = uiputfile({'*.txt;*.log', 'txt / log files'}, ...
    'file to save to', ...
    [get(h.dv, 'string') filesep 'fastECM.log']);
if (fname)
    
    fid = fopen([dname filesep fname], 'w');
    s = get(h.ov, 'string');
    
    if ~iscell(s)
        s = {s};
    end % if
    
    for i = 1:length(s)
        fprintf(fid, '%s\n', s{i});
    end % for
    
    fclose (fid);
    
end % if fname

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% clear the selection of files
%

function clrfiles(h);

set(h.fv, 'string', '', 'value', 0);
clr(h); % also clear log

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% clear the output log
%

function clr(h);

set(  h.ov, 'string', '<empty>', ...
    'value', 1);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% clear the selected mask file
%

function clrmask(h);

set(h.mv, 'string', ' no mask', 'value', 0);
set(h.md, 'string', '' );

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% clear the selected atlas file
%

function clratlas(h);

set(h.av, 'string', ' no atlas', 'value', 0);
set(h.ad, 'string', '' );

return
