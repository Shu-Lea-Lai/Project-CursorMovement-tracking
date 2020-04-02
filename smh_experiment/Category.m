classdef Category
  %   A stimulus category with properties that define the construction of
  %   stimuli within the experiment
  properties
    % what to call it
    name = {};
    % stimulus properties
    
    props = struct();
    cResponse = {};
    
    timingMode = '';
    
    % % % stimulus info % % %
    stimulusInfo = struct();
  end
  methods
    function C = Category(catNum)
      % % % stimulus info % % %
      C.stimulusInfo = C.conf_stimulus();      
      
      % % % category properties % % %
      switch catNum
        case 1
          % stimulus properties, category (A)
          C.props.propShadeMean = [.25];%,.75];
          C.props.propShadeStd = [.1];%,.1];
          C.props.propColorMean = [.25];%,.75];
          C.props.propColorStd = [.1];%,.1];
          C.cResponse = 'F';
        case 2
          % stimulus properties, category (B)
          C.props.propShadeMean = [.25];%,.75];
          C.props.propShadeStd = [.1];%,.1];
          C.props.propColorMean = [.75];%,.25];
          C.props.propColorStd = [.1];%,.1];
          C.cResponse = 'J';
      end
      % ensure that some, but not all pixels are shaded
      C.props.minShaded = 1;
      C.props.maxShaded = (C.stimulusInfo.gridSize).^ 2 - 1;
      C.props.minColored = 0;
      C.props.maxColored = (C.stimulusInfo.gridSize).^ 2 - 1;
      
    end
    
    
  end
  
  methods(Static)    
    %% STIMULUS
    function info = conf_stimulus()
      info.typeName = 'pixelGrid';
      info.gridSize = 10;
      info.Colors = [... % from https://uigradients.com/#Politics
        hex2rgb('f44336');  ... %red
        hex2rgb('2196f3')];     %blue
    end
  end
end
