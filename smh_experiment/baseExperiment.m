  % counterbalance object identities

classdef baseExperiment < handle
  properties
    DESIGN = Design.empty()
    
    DEVICES = struct();
    
    % % % INFO % % %
    info = struct();
  end
  methods
    function E = baseExperiment
      %% Counter-Balance;
      E.info.nOrders = 2;
      E.info.orderSeeds = [111,222];
      
      desiredDir = mfilename('fullpath');
      desiredDir = strrep(desiredDir,[mfilename filesep mfilename],[mfilename filesep]);
      E.info.runLocation = desiredDir;
    end
    
    function E = init_Design(E,S)
      % set seed
      rng(S.config.seed);
        
      cd(E.info.runLocation)
      
      %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%% DESIGN %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Run the experiment      
      % % % numbers % % %
      E.DESIGN = Design();
      E.DESIGN.N.conditions = 1;            
      E.DESIGN.N.blocks = 1;
      E.DESIGN.N.trials = 10;
      
      
      %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%% CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Define conditions
      C{1}.conditionName = 'test';      
      C{1}.Phase(1).stimulus = Stimulus('text');
      C{1}.Phase(1).duration = 1;
      C{1}.Phase(1).device = [];
      C{1}.Phase(2).stimulus = Stimulus('text','test');
      C{1}.Phase(2).duration = 1;
      C{1}.Phase(2).device = [];      
      
      
       
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Setup
      % % % instructions % % %
      E.DESIGN.INSTRUCTIONS = E.init_Instructions(S); % debug
      
      % % % make conditions % % %
      E.DESIGN.CONDITIONS = E.make_Conditions(C);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%% DEVICES %%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Devices
      E.DEVICES.KEYBOARD = Device_Input_Keyboard();
      E.DEVICES.MOUSE = Device_Input_Mouse();
      
       
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % for each block
      for bb = 1:E.DESIGN.N.blocks
        % generate a block
        E.DESIGN.BLOCK(bb) = E.blockMake(E.DESIGN,bb);
        
        % create each trial within the block
        for tt = 1:E.DESIGN.BLOCK(bb).N.trials
          
          % get the design info for this trial
          trialInfo = E.DESIGN.CONDITIONS(E.DESIGN.BLOCK(bb).blockInfo.trialList(tt)).settings;
          
          % create the trial
          E.DESIGN.BLOCK(bb).TRIALS(tt) = E.trialMake(E.DESIGN.designInfo,trialInfo,tt);
        end
      end
      
      
    end
    
    function S = run(E,varargin)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%% CHECK DIR %%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      cd(E.info.runLocation)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%% SUBJECT %%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Import (if provided)
      if nargin == 2
        S = varargin{1};
      else
        E = Experiment('baseExperiment');
        S = E.subs_add;
        % run the subject
        disp(['Running Subject: ' num2str(S.subjectNumber) ' in Order: ' num2str(S.config.order)])
        E.subs_run(S);
        return
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % % % open a window if one doesn't exist % % %
      Device_Output_Display.close_all;
      E.DEVICES.WINDOW = Device_Output_Display(Experiment.init_device_window('full'));
      
      
      Experiment.devices_open(S,E.DEVICES);
      
      global globalOkay
      Experiment.global_escape;
      E.DESIGN.run_Design(S,E.DEVICES);
      S.META.earlyExit = ~globalOkay;
      
    end
    
    
  end
  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function B = blockMake(D,blockNum)
      % create a Block object
      B = Block(blockNum);
      
      % tell it how many trials it should contain
      B.N.trials = D.N.trials;
      
      % how many conditions are there in the experiment?
      nC = length(D.CONDITIONS);
      
      % generate a random order      
      trialOrder = randi(nC,D.N.trials,1);
      
      % define a trial list
      B.blockInfo.trialList = trialOrder;
      
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function T = trialMake(designInfo,trialInfo,trialNum)
      % initialize a Trial object
      T = Trial(trialNum);
      
      % generate the phases
      for pp = 1:length(trialInfo.Phase)
        T.PHASES(pp) = Phase(...
          trialInfo.Phase(pp).stimulus,...
          trialInfo.Phase(pp).duration,...
          trialInfo.Phase(pp).device);
      end
      
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function C = make_Conditions(cList)
      for cc = 1:length(cList)
        C(cc) = Condition(cList{cc});
      end
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% INSTRUCTIONS %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function I = init_Instructions(S)
      I = [];
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% GETS %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function C = get_combined(D)
      C.X = cat(2,D.mouseX{:});
      C.Y = cat(2,D.mouseY{:});
      C.T = cat(2,D.mouseT{:});
      C.B = cat(2,D.mouseB{:});
    end
    
   
    
  end
  
  
  
  
end