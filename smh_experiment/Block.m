classdef Block < handle
  properties
    % % % Number % % %
    blockNum = [];
    
    % % % N % % %
    N = struct('trials',0);
    
    % % % outputs % % %
    outputs = struct();
    
    % % % phase (intro screen) % % %
    PHASES = Phase.empty();
    
    % % % trials container % % %
    TRIALS = Trial.empty();
    
    % % % blockInfo % % %
    blockInfo = struct('outputTable',table.empty);
  end
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% CREATION %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make BLock
    function B = Block(blockNum)
      B.blockNum = blockNum;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run the block
    function B = run_Block(B,S,varargin)
      %% ESCAPE SEQUENCE
      global globalOkay
      
      %% OPTIONAL ARGUMENTS
      % % % DEVICES % % %
      % if the user supplies a list of devices (e.g. PTB Window or
      % keyboard), use those. Otherwise, create a WINDOW at least, for use.
      DEVICES = [];
      if nargin > 1
        DEVICES = varargin{1};
      else
        DEVICES.WINDOW = Device_Output_Display;
      end
      
      
      
      
      %% RUN INTRO
      % % % INTRO PHASES % % %
      % if the user provided some phases to be run at the start of each
      % block, run those first.
      for pp = 1:length(B.PHASES)
        P = B.PHASES(pp);
        
        B.PHASES(pp) = P.run_Phase(DEVICES);
        
        if ~globalOkay
          return
        end
      end
      
      
      %% RUN TRIALS
      % % % TRIALS % % %
      % then run all of the trials for this block.
      B.N.trials = length(B.TRIALS);
      for tt = 1:B.N.trials
        S.DATA.meta.trialNumber = tt;
        
        T = B.TRIALS(tt);
        
        % clear buffers of inputs before starting
        B.device_reset(DEVICES);
        B.TRIALS(tt) = T.run_Trial(S,DEVICES);
        
        B.storeTrial(T,S);
        
        if ~globalOkay
          return
        end
        
      end
      
    end
  end
  methods (Static)
    function device_reset(D)
      dNames = fieldnames(D);
      for dd = 1:length(dNames)
        if strcmp(Experiment.device_type(D.(dNames{dd})),'mouse')
          D.(dNames{dd}).data_reset
        end
      end
    end
    
    function storeTrial(T,S)
      
      % % % META % % %
      data_meta = checkTable(S.DATA.meta);
      
      % % % TRIAL/CONDITION INFO % % %
      data_trial = checkTable(T.trialInfo.outputTable);
      Experiment.store_data([data_meta,data_trial],S,'trial');
      
      % % % PHASE % % %
      data_phase = S.DATA.phase.table;
      phaseOut = [data_phase];
      Experiment.store_data(phaseOut,S,'phase');
      
      
      % % % KEYBOARD % % %
      data_keyboard = S.DATA.keyboard.table;
      keyOut = [data_keyboard];
      Experiment.store_data(keyOut,S,'keyboard');
      
      
      % % % MOUSE % % %
      if S.DATA.use.mouse
        data_mouse = S.DATA.mouse.table;
        mouseOut = [data_mouse];
        Experiment.store_data(mouseOut,S,'mouse')
      end
      
      % % % write the data % % %
      %       S.DATA.mouse = Experiment.store_data(mouseOut,S.DATA.mouse);
      
      
      % % % table check % % %
      % given a structure or table, makes sure it's a table
      function T = checkTable(I)
        if ~isempty(I)
          switch class(I)
            case 'struct'
              T = struct2table(I);
            case 'table'
              T = I;
          end
        else
          T = table.empty();
        end
      end
      
      
      
      
    end
    
    
    
    
  end
end
