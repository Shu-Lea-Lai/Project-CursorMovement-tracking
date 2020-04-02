classdef Phase < handle
  properties
    phaseNum = 1;
    phaseName = '';
    
    % % % STIMULI % % %
    STIMULI = Stimulus.empty
    
    % % % TIMING % % %
    TIMING = struct(...
      'onsetTime',nan,...
      'offsetTime',nan,...
      'earlyTermination',0,...
      'duration', struct('desired',1000,'actual',nan,'error',0));
    
    % % % collect responses % % %
    INPUTS = struct(...
      'waitForResponse',0,...
      'device','',...
      'type','',...
      'allowed',[]);
    
    EVENTS = struct('keyInfo',struct(...
      'keyIsDown', false,...
      'secs', nan,...
      'keyCode',nan,...
      'deltaSecs',nan));
  end
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% CREATE %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make a new phase
    function P = Phase(S,D,varargin)
      for ss = 1:length(S)
        P.STIMULI(ss) = S(ss);
        P.TIMING.duration.desired = D;
      end
      
      % % % add device input as alternative exit condition % % %
      if nargin > 2
        P.INPUTS.waitForResponse = 1;
        
        for ii = 1:size(varargin{1},1)
          P.INPUTS(ii).device = upper(varargin{1}(ii,1));
          P.INPUTS(ii).type = varargin{1}(ii,2);
          P.INPUTS(ii).allowed = varargin{1}(ii,3);
        end
       
      end
     
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run the Phase
    function P = run_Phase(P,DEVICES)
      
      % % % draw it % % %
      for ss = 1:length(P.STIMULI)
        P.STIMULI(ss).run_Stimulus(DEVICES);
      end
      % % % show it % % %
      P.TIMING.onsetTime = DEVICES.WINDOW.Flip;
      % % % wait for it % % %
      P.wait_toContinue(DEVICES);
      % % % record when it's done % % %
      P.TIMING.offsetTime = GetSecs;
      P.TIMING.duration.actual = (P.TIMING.offsetTime - P.TIMING.onsetTime);
      P.TIMING.duration.error = P.TIMING.duration.actual - P.TIMING.duration.desired;
      
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% WAITING %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Wait to continue
    function wait_toContinue(P,DEV)
      %% ESCAPE SEQUENCE
      global globalOkay
      GG = Experiment.global_escape; % global escape seq
      
      %% INITIALIZE CHECKS
      check.keys = 27;
      check.mouse.type = '';
      check.mouse.boxes = [];
      check.eyes = [];
      
      %% INITIALIZE
      time = P.TIMING.duration.desired;
      useKeyboard = 0; keyOKAY = 1; keyInfo = struct('keyIsDown', false,...
      'secs', nan,...
      'keyCode',nan,...
      'deltaSecs',nan);
      useMouse = 0; mouseOKAY = 1; mouseInfo = struct();
      useEyetracker = 0; eyeOKAY = 1; eyeInfo = struct();
      
      
      
      
      
      % % % input devices % % %
      devNames = fieldnames(DEV);
      for ii = 1:length(P.INPUTS)
        if ~isempty(P.INPUTS(ii).device)
          % what time of device is allowed to respond?
          responseDevice = P.INPUTS(ii).device;
          
          thisDevice = DEV.(responseDevice{:});
          % what class is it?
          deviceClass = class(thisDevice);
          
          % % % configure based on class % % %
          switch deviceClass
            case 'Device_Input_Keyboard'
              % turn it on
              useKeyboard = 1;
              
              % convert allowed keys from cell to double if needed
              if iscell(P.INPUTS(ii).allowed)
                AK = cell2mat(P.INPUTS(ii).allowed);
              else
                AK = P.INPUTS(ii).allowed;
              end
              check.keys = [check.keys,AK];
              
              
              K = thisDevice;
              
              
              % % % mouse % % %
            case 'Device_Input_Mouse'
              useMouse = 1;
              M = thisDevice;
              
              check.mouse.type = [check.mouse.type,P.INPUTS(ii).type];
              check.mouse.boxes = [check.mouse.boxes;P.INPUTS(ii).allowed];
              
              % % % eye tracker % % %
            case 'Device_Input_Eyetracker'
              useEyetracker = 1;
              ET = thisDevice;
          end
        end
      end
      
      % % % turn on the keyboard no matter what so we can escape % % %
      if ~useKeyboard
        K = Device_Input_Keyboard;
        useKeyboard = 1;
      end
      
      % % % check timing % % %
      % if desired time is nan, then only move forward with keypress
      if isnan(time)
        time = 1e7; % excessively large length of time
      end
      
      % % % LOOP % % %
      startTime = GetSecs;
      while [(GetSecs - startTime) <= time] && keyOKAY && mouseOKAY && eyeOKAY
        % % % check keyboard % % %
        if useKeyboard
          % % % user inputs % % %
          keyInfo = K.get_Keys;
          keyOKAY = ~K.compareKeys(keyInfo.keyCode,check.keys);
          
          
          % % % global escape sequence % % %
          globalOkay = ~K.compareKeys(keyInfo.keyCode,GG);
          if ~globalOkay
            P.TIMING.earlyTermination = 1;
            storeEvents(P,keyInfo,mouseInfo,eyeInfo)
            return
          end
        end
        % % % check mouse % % %
        if useMouse
          % poll the mouse
          M.get_Mouse;
          
          % check mouse against hit/in
          wasClicked = 0;
          wasIn = 0;
          for cc = 1:length(check.mouse.type)
            switch check.mouse.type{cc}
              case 'click'
                wasClicked = wasClicked + M.clickTest(...
                  check.mouse.boxes{cc,:});
                [M.mouseInfo.x(M.i-1),M.mouseInfo.y(M.i-1)];
                
              case 'over'
                wasIn = wasIn + M.hitTest(...
                  check.mouse.boxes{cc,:});
                [M.mouseInfo.x(M.i-1),M.mouseInfo.y(M.i-1)];
                check.mouse.boxes{cc,:};
            end
          end
          if wasIn || wasClicked
            mouseOKAY = 0;
          end
        end
        % % % check eyetracker % % %
        if useEyetracker
        end
      end
      
      storeEvents(P,keyInfo,mouseInfo,eyeInfo)
      
      function storeEvents(P,keyInfo,mouseInfo,eyeInfo)
        P.EVENTS.keyInfo = keyInfo;
        P.EVENTS.mouseInfo = mouseInfo;
        P.EVENTS.eyeInfo = eyeInfo;
      end
    end
  end
  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% MISC FCNS %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Seconds to MS
    function t = s2ms(t,varargin)
      % % % convert times in seconds to milliseconds % % %
      % can force past smart checking
      force = 0; if nargin > 1 && strcmp(varargin{1},'f'); force = 1; end
      if t <= 10 || force
        t = t .* 1000;
      else
        disp(['TIME: ' num2str(t) ' NOT CONVERTED TO MS'])
      end
    end
    %% MS to Seconds
    function t = ms2s(t,varargin)
      % % % convert times in milliseconds to seconds % % %
      % can force past smart checking
      force = 0; if nargin > 1 && strcmp(varargin{1},'f'); force = 1; end
      if t >= 10 || force
        t = t ./ 1000;
      else
        disp(['TIME: ' num2str(t) ' NOT CONVERTED TO SECONDS'])
      end
    end
    
    
  end
end