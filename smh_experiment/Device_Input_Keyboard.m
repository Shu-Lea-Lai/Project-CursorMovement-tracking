classdef Device_Input_Keyboard < handle
  properties
    % % % handle % % %
    kbNum = 0;
  end
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% SELF %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Create
    function K = Device_Input_Keyboard()      
    end
    
    %% Collect Keypresses
    function keyInfo = get_Keys(K)
      [keyInfo.keyIsDown, keyInfo.secs, keyInfo.keyCode, keyInfo.deltaSecs] = KbCheck(K.kbNum);
    end    
    
    
    
    
  end
  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% CONVERT %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% KbName
    function keyNum = nameNum(keyName)
      % converts key names in strings to numbers for faster comparison
      if ischar(keyName)
        keyNum = KbName(keyName);
      else
        keyNum = keyName;
      end      
    end
    
    function keyName = numName(keyNum)
      % converts key numbers back to names
      if isnumeric(keyNum)
        keyName = KbName(keyNum);
      else
        keyName = keyNum;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% COMPARE %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compare Keys
    function keyMatch = compareKeys(pressed,queried)
      pressedNum = find(pressed);
      keyMatch = any(any(pressedNum == queried'));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% UTILITIES %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function getKeyNames
      startTime = GetSecs;
      while GetSecs - startTime < 3
        [~,~,keyName] = KbCheck;
        KbName(keyName)
      end
    end
  end
end