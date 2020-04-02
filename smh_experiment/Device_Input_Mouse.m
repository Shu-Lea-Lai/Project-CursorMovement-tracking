classdef Device_Input_Mouse < handle
  properties
    % % % handle % % %
    mouseNum = 0;
    % % % index % % %
    i = 1;
    % % % storage % % %
    mouseInfo = struct(...
      'x',nan(1,10000),...
      'y',nan(1,10000),...
      't',nan(1,10000),...
      'b',nan(3,10000));
  end
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% SELF %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Create
    function M = Device_Input_Mouse()
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% COLLECT %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function get_Mouse(M)
      [M.mouseInfo.x(M.i), M.mouseInfo.y(M.i),buttons] = GetMouse();
      if length(buttons) < 3
        M.mouseInfo.b(:,M.i) = [buttons, nan(1,3-length(buttons))];
      elseif length(buttons) > 3
        M.mouseInfo.b(:,M.i) = buttons(1:3);
      else
        M.mouseInfo.b(:,M.i) = buttons;
      end
      M.mouseInfo.t(M.i) = GetSecs;
      M.i = M.i + 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% PRUNING %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Read data from buffer
    function mouseInfo = data_read(M)
      fn = fieldnames(M.mouseInfo);
      for ff = 1:length(fn)
        V = M.mouseInfo.(fn{ff});
        mouseInfo.(fn{ff}) = V(:,1:M.i-1);
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% RESET %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Reset buffer
    function data_reset(M)
      M.mouseInfo = struct(...
        'x',nan(1,10000),...
        'y',nan(1,10000),...
        't',nan(1,10000),...
        'b',nan(3,10000));
      M.i = 1;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% CONTINGENCY %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function yesNo = hitTest(M,box)
      yesNo = [M.mouseInfo.x(M.i-1) >= box(1)] & [M.mouseInfo.x(M.i-1) <= (box(1) + box(3))] & ...
          [M.mouseInfo.y(M.i-1) >= box(2)] & [M.mouseInfo.y(M.i-1) <= (box(2) + box(4))];        
    end
    
    function yesNo = clickTest(M,box)
      if any(M.mouseInfo.b(:,M.i-1))
        yesNo = [M.mouseInfo.x(M.i-1) >= box(1)] & [M.mouseInfo.x(M.i-1) <= (box(1) + box(3))] & ...
          [M.mouseInfo.y(M.i-1) >= box(2)] & [M.mouseInfo.y(M.i-1) <= (box(2) + box(4))];
      else
        yesNo = 0;
      end
    end
    
    
  end
  methods (Static)
    
    
    % define a static function that does the same as hitTest
    function yesNo = static_hitTest(XY,box)
      yesNo = [XY(:,1) >= box(1)] & [XY(:,1) <= (box(1) + box(3))] & ...
          [XY(:,2) >= box(2)] & [XY(:,2) <= (box(2) + box(4))];
    end
    
    
    
    
  end
  
end