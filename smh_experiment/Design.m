classdef Design < handle
  properties
    % % % how many of things % % %
    N = struct();
    
    % % % conditions % % %
    CONDITIONS = Condition.empty();
    
    % % % outputs % % %
    outputs = struct();
    
    % % % subcontainer for BLOCK % % %
    BLOCK = Block.empty();
    
    % % % subcontainer for instructions % % %
    INSTRUCTIONS = Instructions.empty();
    
    % % % misc % % %
    designInfo = struct();
    
  end
  methods
    function D = Design(D)
    end
    
    function D = add_Block(D)
      currentBlocks = length(D.BLOCK);
      D.BLOCK(currentBlocks+1) = Block(currentBlocks+1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function run_Design(D,S,varargin)
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
      
      %% RUN INSTRUCTIONS (IF PRESENT)
      if ~isempty(D.INSTRUCTIONS)
        D.INSTRUCTIONS.run_Instructions(S,DEVICES);
      end
      
      %% RUN BLOCKS
      % Run each of the blocks associated with this design.
      for bb = 1:D.N.blocks
        % % % store data % % %
        S.DATA.meta.blockNumber = bb;
        
        
        B = D.BLOCK(bb);
        D.BLOCK(bb) = B.run_Block(S,DEVICES);
                
        if ~globalOkay
%           D.exit_proc(S,DEVICES);
          return
        end
      end
%       D.exit_proc(S,DEVICES);
    end
  end
  methods (Static)
    
%     function exit_proc(S,DEVICES)      
%       fclose('all');
%     end
%     
%     function write_data(S)
%       % % % keyboard % % %
%       Experiment.write_data(S,'keyboard');
%       
%       % % % mouse % % %
%       Experiment.write_data(S,'mouse');
%       
%       % % % eye % % %
%     end
    
  end
end