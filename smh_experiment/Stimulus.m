classdef Stimulus < handle
  properties (SetObservable)
    % % % meta % % %
    name = '';
    
    % % % type % % %
    type = '';
    
    % % % properties % % %
    PROPS = struct();
  end
  events
    fileNameChange
  end
  methods
    function S = Stimulus(type,varargin)
      % % % add type % % %
      S.type = type;
      % % % assign properties % % %
      switch type
        case 'blank'
          % by default, color = 'black';
          S.PROPS.color = [0 0 0];
          % % % update with inputs % % %
          S.PROPS = S.addOptional(S.PROPS,varargin);
        case 'text'
          % % % defaults % % %
          % text properties
          S.PROPS.str = 'hello world';
          S.PROPS.color = [0 0 0];
          S.PROPS.fontSize = 16;
          % position
          S.PROPS.positionX = 0.5;
          S.PROPS.positionY = 0.5;
          % % % update with inputs % % %
          S.PROPS = S.addOptional(S.PROPS,varargin);
        case 'fileImage'
          % % % file name-change listener % % %
          addlistener(S, 'fileNameChange', @S.load_newFileName);
          % % % defaults % % %
          % color data for drawing
          % file name
          S.PROPS.fileName = 'office_4.jpg';
          % ptb texture
          S.PROPS.ptbTexture = [];
          % position
          S.PROPS.positionX = .5;
          S.PROPS.positionY = .5;
          % % % update with inputs % % %
          S.PROPS = S.addOptional(S.PROPS,varargin);
        case 'matrixImage'
          % % % defaults % % %
          % color data for drawing
          S.PROPS.matrix = rand(10,10);
          S.PROPS.image = rand(100,100);
          % ptb texture
          S.PROPS.ptbTexture = [];
          % position
          S.PROPS.positionX = .5;
          S.PROPS.positionY = .5;
          % % % update with inputs % % %
          S.PROPS = S.addOptional(S.PROPS,varargin);
        case 'sound'
          S.PROPS = struct(...
            'fileName','',...
            'wavData',[],...
            'freq',[],...
            'nChannels',[]...
            );
          if nargin > 1
            S.PROPS = varargin{1};
          end
        case 'video'
        otherwise
          disp(['WARNING: Could not create Stimulus of type: ' type])
      end
    end
    function set.PROPS(S,property)
      FN = fieldnames(property);
      for ff = 1:length(FN)
        % store new and old values for comparison later
        newVal = property.(FN{ff});
        oldVal = '';
        if isfield(S.PROPS, FN{ff})
          oldVal = S.PROPS.(FN{ff});
        end
        S.PROPS.(FN{ff}) = property.(FN{ff});
        % check if filename has changed
        if strcmp(FN{ff},'fileName') && ~strcmp(newVal,oldVal)
          notify(S,'fileNameChange')
        end
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% LOADING %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Load New Input
    function S = load_newFileName(S,varargin)
      switch S.type
        case 'fileImage'
          S = S.load_fileImage();
        case 'sound'
          S = S.load_fileSound();
      end
    end
    
    %% Load fileImage
    function S = load_fileImage(S,varargin)
      S.PROPS.matrix = imread(S.PROPS.fileName);
      S.PROPS.image = S.PROPS.matrix;
    end
    
    %% load fileSound
    function S = load_fileSound(S,varargin)
      [S.PROPS.wavData, S.PROPS.freq] = audioread(S.PROPS.fileName);
      S.PROPS.nChannels = size(S.PROPS.wavData,2);
      try
        PsychPortAudio('Close');
      catch
        disp(['didnt load new ' S.PROPS.filename]);
      end
      S.PROPS.h = PsychPortAudio('Open', [], [], 0, S.PROPS.freq, S.PROPS.nChannels);
      PsychPortAudio('FillBuffer', S.PROPS.h, S.PROPS.wavData');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% RUNNING %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function S = run_Stimulus(S,DEV)
      switch S.type
        case {'text','fileImage','matrixImage'}
          S.draw_Stimulus(DEV.WINDOW);
        case 'sound'
          S.play_Stimulus(DEV.SOUND);
      end
    end
    
    %% Draw a Visual Stimulus
    function draw_Stimulus(S,WIN)
      WIN.draw_Stimulus(S);
    end
    
    %% Play an Auditory Stimulus
    function play_Stimulus(S,SOUND)
      SOUND.play_Sound(S);
    end
  end
  methods (Static)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Use varargins
    function P = addOptional(P,ins)
      % if additional arguments are given during creation, mofidy the
      % defaults
      fn = fieldnames(P);
      for i = 1:length(ins)
        % allow user to leave some as deffaults by making the inputs empty
        if ~isempty(ins{i})
          P.(fn{i}) = ins{i};
        end
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% TOOLS %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    
    
  end
end