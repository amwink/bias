function P = spm_get(varargin)

% This programm converts the spm_get command into a spm_select
% command, which is the file-selector in SPM5. It allows old scripts
% and functions to be used within spm5 and is backward compatible.
%
% It transforms spm<5  syntax:    spm_get(n,filt,Prompt,NewWDir,CmdLine);
% into          spm5   syntax:    spm_select(n,typ,mesg,sel,wd,filt,frames);
%
% Usage:
%========
% P = spm_get(n,filt,Prompt,NewDir,CmdLine)
% 
% n         = number of files to be selected (neg. numbers for directories)
% filt      = filename filter
% Prompt    = Message for the user
% NewDir    = directory to start looking in
% CmdLine   = whether to use GUI (1) or command line (0); unused in spm5
%
% All these input options are optional
% 
% The program converts it into:
%
% P = spm_select(n,type,Prompt,Sel,wd,filt,frames);
%
% with:
%
% n      =  number of files to be selected (only pos.)
%           (defaults to 'Inf' if not specified)
% type   =  either 'image', 'dir', or the converted filt
%           (defaults to 'any' if not specified)
% Prompt =  contains the original Prompt with the determined type in brackets 
%           (defaults to 'Select file' if not specified)
% Sel    =  list of already selected files (not available in spm<5, so
%           defaults to [] in any case)
% wd     =  directory to start in, from NewDir
%           (defaults to pwd if not specified)
% filt   =  prefix to look for, derived from the filt input
%           (defaults to '.*' if not specified)
% frames =  image frame numbers to include (has no equivalent in spm<5, so
%           defaults to 1 in any case)
%_______________________________________________________________________
% Copyright (C) 2006 Karsten Specht, IBMP, Bergen, Norway
%                    with input from Marko Wilke, Tuebingen, Germany

%
% $Id: spm_get.m 1.3 2006-09-26  Karsten & Marko $

% get, check version (conversion may not be necessary at all)
%--------------------------------
  ver = spm('ver');
  if strcmp(ver,'SPM2') == 1 | strcmp(ver,'SPM99') == 1
	  fprintf(['   ' ver ' detected, reverting to old syntax ...' '\n']);


	% find the original spm_get (likely the second one in the path) and go there
	%----------------
	%  W = which('spm_get.m','-all');
	%  [pth nm e] = fileparts(W{2,:});
	%  a = pwd;
	%  cd(pth);

    % Go to the spm-folder. spm_get should be there
    %---------------------------------------------
    a = pwd;
    cd(spm('dir'));


	% convert input to output
	%----------------
	  if nargin == 5
		n = varargin{1};     filt = varargin{2};     Prompt = varargin{3};           NewDir = varargin{4};     CmdLine = varargin{5};
	  elseif nargin == 4
		n = varargin{1};     filt = varargin{2};     Prompt = varargin{3};           NewDir = varargin{4};     CmdLine = 0;
	  elseif nargin == 3
		n = varargin{1};     filt = varargin{2};     Prompt = varargin{3};           NewDir = a;               CmdLine = 0;
	  elseif nargin == 2 
		n = varargin{1};     filt = varargin{2};     Prompt = 'Select files...';     NewDir = a;               CmdLine = 0;
	  elseif nargin == 1
		n = varargin{1};     filt = '*.*';           Prompt = 'Select files...';     NewDir = a;               CmdLine = 0;
	  elseif nargin == 0
		n = Inf;             filt = '*.*';           Prompt = 'Select files...';     NewDir = a;               CmdLine = 0;
	  end;


	% run original spm_get
	%----------------
	  P = spm_get(n,filt,Prompt,NewDir,CmdLine);
	  cd(a);
	  return
      
  else

	% Set the default values
	%--------------------------------
	  sel    = [];
	  frames = '1';

	
	% Determine the number of files
	%--------------------------------
	  if nargin > 0
		n = varargin{1};
	  else
		n = Inf;
	  end;

	
	% Determine the type of selection and the filter
	%---------------------------------
	  if nargin > 1
		filt = varargin{2};
	
 		if ~isempty(strfind(filt,'*.mat'))
 	        	type = 'mat';
                filt = spm_str_manip(filt,'s');
        elseif ~isempty(strfind(filt,'*.img')) | ~isempty(strfind(filt,'*.nii')) | ~isempty(strfind(filt,'*.mnc'))
 	        	type = 'image';
                filt = spm_str_manip(filt,'s');
        elseif ~isempty(strfind(filt,'image')) | ~isempty(strfind(filt,'IMAGE'))
 	        	type = 'image';
           		filt = '.*';
        else
                filt = spm_str_manip(filt,'s');
 			    type = 'any';
        end;

      else
		type = 'any';
		filt = '.*';
	  end;
	
	
	% Determine the Prompt
	%-------------------------
	  if nargin > 2
		Prompt = varargin{3};
		Prompt = [Prompt ' (' type ')'];
	  else
		Prompt = 'Select file';
	  end;

	
	% Special case: Directory
	%-------------------------
	  if n < 0
		n      = abs(n);
		type   = 'dir';
		Prompt = 'Select directory';
	  end;


	% Get directory
	%--------------------------------
	  if nargin > 3
		wd = varargin{4};
	  else
		wd = pwd;
	  end;
	
	% Run spm_select
	%----------------
	  P = spm_select(n,type,Prompt,sel,wd,filt,frames);
end;