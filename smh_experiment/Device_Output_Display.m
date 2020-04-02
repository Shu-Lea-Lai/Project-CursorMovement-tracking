classdef Device_Output_Display < handle
  properties
    % % % handle % % %
    H = [];
    % % % rectangle % % %
    rect = [];
    % % % properties % % %
    props = struct();
  end
  methods
    function W = Device_Output_Display(properties)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%% OPEN  %%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Screen('Preference','SkipSyncTests',properties.skipSync);
      % % % open the screen % % %
      % if the skipSync flag is set to 0, we really care about timing.
      % % Unfortunately, opening a window with this setting can often fail
      % % the first time, but succeed if you try again. As a result, we're
      % % going to try more than once, if it fails.
      W.props.tries = 0;
      W.H = [];
      while isempty(W.H) && W.props.tries < 3
        try
          [W.H,W.rect] = Screen('OpenWindow',...
            properties.screenNum,... % SCREEN NUM
            properties.screenColor,... % COLOR
            properties.screenRect); %  RECT
        end
        W.props.tries = W.props.tries + 1;
      end
      
      % if the screen refuses to open anyway, turn on Sync and try again
      Screen('Preference','SkipSyncTests',1);
      [W.H,W.rect] = Screen('OpenWindow',...
        properties.screenNum,... % SCREEN NUM
        properties.screenColor,... % COLOR
        properties.screenRect); %  RECT
      W.props = properties;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% DISPLAYING %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generic Draw Function
    function S = draw_Stimulus(W,S)
      switch S.type
        case 'text'
          W.draw_Text(S.PROPS);
        case {'fileImage','matrixImage'}
          S.PROPS.ptbTexture = W.create_Texture(S.PROPS.image);
          S.PROPS.rect = W.compute_rect(S);
          W.draw_Texture(S.PROPS);
      end
    end
    
    %% DrawText
    function draw_Text(W,text)
      % % % configure text size % % %
      setTextSize(W,text);
      % % % get coordinates % % %
      extent = textExtent(W,text);
      % % % where to begin the text % % %
      P = getObjectLocations(W,text,extent);
      Screen('DrawText', ...
        W.H,... % handle
        text.str,... % string
        P(1),P(2),... % position
        text.color);
    end
    
    %% Create a Texture (image)
    function T = create_Texture(W,matrix)
      T = Screen('MakeTexture',W.H,matrix);
    end
    
    %% Draw a texture
    function draw_Texture(W,PROPS)
      Screen('DrawTexture',W.H,PROPS.ptbTexture,[],PROPS.rect);
    end
    
    %% Flip
    function flipTime = Flip(W)
      flipTime = Screen('Flip',W.H);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% IMAGE SIZING %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Rect for Images
    function R = compute_rect(W,S)
      xCntrPix = round(S.PROPS.positionX * W.rect(3));
      yCntrPix = round(S.PROPS.positionY * W.rect(4));
      
      imSize = size(S.PROPS.image);
      
      R = [xCntrPix - (imSize(2)/2), yCntrPix - (imSize(1)/2), ...
        xCntrPix + (imSize(2)/2), yCntrPix + (imSize(1)/2)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% TEXT SIZING %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Text Extent
    function E = textExtent(W,text)
      E = Screen('TextBounds', W.H, text.str, text.positionX, text.positionY);
    end
    %% Scale Text Position
    function text = scaleTextPosition(W,text)
      % if text positions are 0-1, scale to the size of the window
      if text.positionX < 1
        text.positionX = text.positionX * W.props.screenRect(3);
      end
      if text.positionY < 1
        text.positionY = text.positionY * W.props.screenRect(4);
      end
    end
    %% Text Size
    function setTextSize(W,text)
      % % % Text Size % % %
      Screen('TextSize',W.H,text.fontSize);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% LOCATIONS %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function P = getObjectLocations(W,obj,E)
      P(1) = (W.rect(3) * obj.positionX) - (E(3)/2);
      P(2) = (W.rect(4) * obj.positionY) - (E(4)/2);
    end
    
    
    
    
  end
  
  methods(Static)
    %% Close all
    function close_all
      % % % close all open windows % % %
      sca
    end
    
    %% Open a test window
    function info = open_demo(varargin)
      % which screen to use?
      if nargin == 0
        screenNum = 0;
      else
        screenNum = varargin{1};
      end
      
      % skip sync?
      info.skipSync = 1;
      % PTB screen configuration
      info.groot = groot;
      % which screen to use
      info.screenNum = screenNum;
      % what color to fill screen with?
      info.screenColor = [255 255 255];
      % screen size/location
      info.screenRect = [10 10 750 750];
    end
    
    
    %% Open a full window
    function info = open_full(varargin)
      % which screen to use?
      if nargin == 0
        screenNum = 0;
      else
        screenNum = varargin{1};
      end
      
      % skip sync?
      info.skipSync = 0;
      % PTB screen configuration
      info.groot = groot;
      % which screen to use
      info.screenNum = screenNum;
      % what color to fill screen with?
      info.screenColor = [255 255 255];
      % screen size/location
      info.screenRect = [0 0 info.groot.ScreenSize(3) info.groot.ScreenSize(4)];
    end
  end
end