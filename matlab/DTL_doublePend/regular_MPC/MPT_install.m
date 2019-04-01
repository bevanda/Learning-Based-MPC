%% installation of tbxmanager with all submodels required for MPT
%

clc;
disp('----------------------------------------------');
disp('Installation of MPT using the Toolbox manager.');
disp('----------------------------------------------');
disp(' ');
fprintf(['Choose the installation directory where to install the Toolbox manager.\n',...
      'A new folder "tbxmanager" is going to be created in the specified location.\n',...
      'If you do not specify the folder, the Toolbox manager will be installed in the current directory.\n']);

% get the installation folder
default_dir = pwd;
c = uigetdir(pwd);
if isequal(c,0);y
    fprintf(['No directory has been provided.\n',... 
        'Installing the toolbox manager in the current directory "%s"?\n'],default_dir);
    c = default_dir;
end
 
% create a new directory in that folder
d = [c,filesep,'tbxmanager'];
if isequal(exist(d,'dir'),7)
    error('The installation directory "%s" already exists.\nPlease, remove or rename the folder or change the installation path.',d);
end
disp('Creating the directory "tbxmanager".');
out = mkdir(d);
if ~out
    error(['An error appear when trying to create the folder "%s".\n',...
          'Please, install the Toolbox manager manually.'],c); 
end

% enter that directory
cd(d);

% remove MPT2 or YALMIP
disp(' ');
disp('Removing toolboxes that may conflict with MPT from the Matlab path.');
rmpath(genpath(fileparts(which('mpt_init'))));
rmpath(genpath(fileparts(which('yalmipdemo'))));


% download the tbxmanager
disp(' ');
disp('Downloading the Toolbox manager from the internet.');
[f, c] = urlwrite('http://www.tbxmanager.com/tbxmanager.m', 'tbxmanager.m');
rehash;

if isequal(c,0)
    error('Could not download the Toolbox manager from the internet. The installation cannot continue.');
end

% install all required modules
tbxmanager install mpt mptdoc cddmex fourier glpkmex hysdel lcp sedumi yalmip 

% create the initialization file to set the path 
disp(' ');
disp('Creating the initialization file "startup.m".');
p = which('startup.m');
if isempty(p)
    p = [d,filesep,'startup.m'];
end
fid = fopen(p,'a');
if isequal(fid,-1)
    error(['Could not modify the initialization file "startup.m".',...
           'Edit this file in the folder "%s" manually and insert there the line:  tbxmanager restorepath.'],p);
end
fprintf(fid,'tbxmanager restorepath\n');
fclose(fid);
disp('File has been created.');

% get back to the original directory
cd(default_dir);

% add path to tbxmanager
disp(' ');
disp('Adding path to Matlab.');
addpath(d);

% save path for future
disp(' ');
disp('Saving path for future sessions.');
status = savepath;

if status
    fprintf('Could not save the path to a default location,\nplease provide a location where you want to save the path.');
    cn = uigetdir(pwd);
    if isequal(cn,0)
        disp(' ');
        fprintf('No directory specified, saving the path to the current directory "%s".\n\n',default_dir);
        cn = default_dir;
    end
    sn = savepath([cn,filesep,'pathdef.m']);
    if sn
        error(['Could not save the path automatically.\n',...
            'Please, open the "Set Path" button in the Matlab menu and save the path manually to some location.']);
    end
end

disp(' ');
disp('Installation finished.');
disp('Next time you start Matlab the toolboxes will be automatically initialized.');

% initialize MPT
disp(' ');
disp('Initializing the MPT.')
mpt_init;
