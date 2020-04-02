classdef Device_Output_Sound < handle
  properties
    % % % handle % % %
    h = [];
    
    % % % BUFFERS % % %
    BUFFERS = Stimulus.empty()
    
    
  end
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% OPEN  %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialize Sound
    function S = Device_Output_Sound()
      InitializePsychSound;
      try
        S.h = PsychPortAudio('Open');
      catch
        % sound may still be turned on, causing an error, so disable and
        % tray again
        PsychPortAudio('Close')
        S.h = PsychPortAudio('Open');
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% BUFFERS %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Create a new Buffer
    function add_NewBuffer(S,fileName)
      nBuffers = length(S.BUFFERS); tB = nBuffers + 1;
      S.BUFFERS(tB) = S.create_buffer(fileName);
    end
    
    %% Add an existing Buffer
    function add_ExistingBuffer(S,B)
      nBuffers = length(S.BUFFERS); tB = nBuffers + 1;
      S.BUFFERS(tB) = B;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PLAY %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Play Sound
    function play_Sound(S,B)
      PsychPortAudio('FillBuffer', S.h, B.PROPS.wavData');
      startTime = PsychPortAudio('Start',S.h,1,0,1);

      s = PsychPortAudio('GetStatus', S.h);
      s.Active = 1;
      while s.Active == 1
        s = PsychPortAudio('GetStatus', S.h);
      end
      stopTime = PsychPortAudio('Stop',S.h);
    end
  end
  
  methods (Static)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% BUFFER %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function B = create_buffer(fileName)
      props.fileName = fileName;
      [props.wavData, props.freq] = psychwavread(fileName);
      props.nChannels = size(props.wavData,2);
      
      B = Stimulus('sound',props);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% CLOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Close all Sounds
    function close_all
      PsychPortAudio('Close');
    end
  end
end