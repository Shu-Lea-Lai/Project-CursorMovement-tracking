classdef Instructions < handle
  properties
    % % % PHASES container % % %
    PHASES = Phase.empty;
  end
  
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run the block
    function B = run_Instructions(I,S,varargin)
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
      
      
      %% KEYBOARD
      % We will need to control stimulus presentation somehow, and we will
      % use a keyboard by default. If one isn't enabled, create one
      if ~isfield(DEVICES,'KEYBOARD')
        DEVICES.KEYBOARD = Device_Input_Keyboard;
      end
      
      
      %% RUN INTRO
      % % % INTRO PHASES % % %
      % if the user provided some phases to be run at the start of each
      % block, run those first.
      pp = 1;
      while pp <= length(I.PHASES)
        P = I.PHASES(pp);
        
        I.PHASES(pp) = P.run_Phase(DEVICES);
        
        
        % % % wait for a key-press % % %
        pp = Instructions.phase_next(P,pp);
        KbReleaseWait
        
        if ~globalOkay
          return
        end
        
        
        
      end
      
    end
    
  end
  
  methods (Static)
    
    function nextP = phase_next(P,pp)
      % set default value in case of early termination
      nextP = pp;
      % gather key events
      keyEvents = P.EVENTS.keyInfo;
      keyName = KbName(keyEvents.keyCode);
      
      % some slides advance automatically
      if isempty(keyName)
        nextP = pp + 1;
        % otherwise, see if we want to go forward or backwards
      else
        switch keyName
          case {'space',32,'right',39}
            nextP = pp + 1;
          case {'left',37}
            nextP = max(1,pp - 1);
        end
      end
    end
  end
end