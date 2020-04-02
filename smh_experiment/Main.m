% MAIN
% function to run experiments from
sca
close all

% % % GLOBAL ESCAPE SEQUENCE % % %
global ESC
ESC = 0;
%% Setep / Init

% % % add dependencies % % %
dependDir = [pwd filesep 'depend'];
addpath(dependDir)


%% Create the Experiment
% create Experiment object
E = Experiment();


% % %  % % %
f.h = figure(1);
f.ui(1) = uicontrol();
f.ui(1).Parent = f.h;
f.ui(1).Style = 'pushbutton';
f.ui(1).Callback = {@runButton,E};

%% Run the Experiment
% E.run_subject


function runButton(h,o,E)
global ESC
while ESC == 0
  E.run_subject
end
keyboard
end