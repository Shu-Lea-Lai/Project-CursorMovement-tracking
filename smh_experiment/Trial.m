classdef Trial < handle
  properties
    % % % Number % % %
    trialNum = [];
    
    % % % N % % %
    N = struct('phases',0);
    
    % % % outputs % % %
    outputs = struct();
    
    % % % PHASES container % % %
    PHASES = Phase.empty;
    
    % % % catch-all for other items % % %
    trialInfo = struct('outputTable',table.empty);
  end
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% CREATION %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make Trial
    function T = Trial(trialNum)
      T.trialNum = trialNum;
    end
    %% Add a new Trial
    function T = add_Phase(T)
      currentPHASES = length(T.PHASES);
      T.PHASES(currentPHASES+1) = Trial(currentPHASES+1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run the trial
    function T = run_Trial(T,S,varargin)
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
      
      
      
      %% RUN
      % Run each of the phases associated with this trial.
      T.N.phases = length(T.PHASES);
      for pp = 1:T.N.phases
        % % % wait for a release on all keyboards % % %
        KbReleaseWait;
        
        
        P = T.PHASES(pp);
        P.phaseNum = pp;
        T.PHASES(pp) = P.run_Phase(DEVICES);
        
        if ~globalOkay
          T.storeData(T,DEVICES,S);
          return
        end
      end
      T.storeData(T,DEVICES,S);
      
    end
  end
  methods (Static)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% WRAPPER %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function storeData(T,DEVICES,S)
      for pp = 1:T.N.phases
        P = T.PHASES(pp);
        colHead = T.getPhaseName(P);
        
        S.DATA.meta.trialNumber = T.trialNum;
        % % % phase data % % %
        S.DATA.phase = T.data_basic(P,S.DATA.phase,colHead);
        
        % % % loop through devices % % %
        dFN = fieldnames(DEVICES);
        for dd = 1:length(dFN)
          D = DEVICES.(dFN{dd});
          switch class(D)
            case 'Device_Input_Keyboard'
              keyInfo = P.EVENTS.keyInfo;
              S.DATA.keyboard = T.data_keyboard(keyInfo,S.DATA.keyboard,colHead);
            case 'Device_Input_Mouse'
              if pp == T.N.phases
                S.DATA.mouse = T.data_mouse(D,S.DATA.mouse);
                D.data_reset;
              end
          end
        end
      end
    end
    
    
    %% Phase Name
    function columnHead = getPhaseName(P)
      if ~isempty(P.phaseName)
        columnHead = P.phaseName;
      else
        columnHead = ['P' num2str(P.phaseNum)];
      end
      % append an underscore
      columnHead = [columnHead '_'];
    end
    
    
    %% RECORD TRIAL INFO
    % For each phase, record basic information like on/off-set times and
    % any event data that was collected.
    function DATA = data_basic(P,DATA,colHead)
      
      % % % add basic timing info % % %
      DATA.table.([colHead 'onsetTime']) = P.TIMING.onsetTime;
      DATA.table.([colHead 'offsetTime']) = P.TIMING.offsetTime;
      DATA.table.([colHead 'earlyTermination']) = P.TIMING.earlyTermination;
      DATA.table.([colHead 'duration_desired']) = P.TIMING.duration.desired;
      DATA.table.([colHead 'duration_actual']) = P.TIMING.duration.actual;
      DATA.table.([colHead 'duration_error']) = P.TIMING.duration.error;
      
    end
    
    
    function DATA = data_keyboard(keyInfo,DATA,colHead)
      
      % % % add keyboard info % % %
      DATA.table.([colHead 'keyboard_keyPressed']) = keyInfo.keyIsDown;
      DATA.table.([colHead 'keyboard_keyTime']) = keyInfo.secs;
      % get name of key if pressed
      keyName = find(keyInfo.keyCode);
      
      % do some managing of keys
      if isempty(keyName)
        % % if a key wasn't pressed, set time = nan
        DATA.table.([colHead 'keyboard_keyTime']) = nan;
        keyName = nan;
      elseif length(keyName) > 1
        % % if they pressed 2 keys at once, simulate "backspace," as this
        % % should not conflict with any study. Still record the time as
        % % nan.
        DATA.table.([colHead 'keyboard_keyTime']) = nan;
        keyName = 8;
      end
      DATA.table.([colHead 'keyboard_keyName']) = keyName;
    end
    
    function DATA = data_mouse(M,DATA)
      % % % extract mouse info % % %
      mouseInfo = M.data_read;
      
      DATA.table.mouseT = {mouseInfo.t};
      DATA.table.mouseX = {mouseInfo.x};
      DATA.table.mouseY = {mouseInfo.y};
      DATA.table.mouseB = {mouseInfo.b};
      
      % % % add ET info % % %
      % to be added
    end
    
  end
end
