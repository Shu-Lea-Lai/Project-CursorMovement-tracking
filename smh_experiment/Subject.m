classdef Subject < handle
  properties
    % % % number % % %
    subjectNumber = [];
    
    % % % experiment Name % % %
    experimentName = '';
    
    % % % meta % % %
    META = struct('didRun',0,'didFinish',0,'runTime','','finishTime','','earlyExit',nan,'computername','');
    
    %     % % % exp data % % %
    %     EXP = struct();
    
    % % % data % % %
    DATA = struct(...
      'use',struct(),...
      'meta',struct(),...
      'trial',struct(...
      'table',table.empty,...
      'tableFull',table.empty(),...
      'columns',table.empty()),...
      'phase',struct(...
      'table',table.empty,...
      'tableFull',table.empty(),...
      'columns',table.empty()),...
      'keyboard',struct(),...
      'mouse',struct(),...
      'eye',struct());
    
    % % % analyses % % %
    analyses = struct();
    
    % % % plots % % %
    plotting = struct();
    
    
    % % % config % % %
    config = struct('order',nan,'seed',nan);
  end
  methods
    %% Initialize
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function S = Subject(experimentName)      
      % % % save the name of the experiment % % %
      S.experimentName = experimentName;      
    end
    
    %% Data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % Initialize % % %
    function data_init(S,DEVICES)
      fileName = [S.experimentName '_Subject' num2str(S.subjectNumber)];
      S.DATA.meta.fileHeader = {fileName};
      S.DATA.meta.subjectNumber = S.subjectNumber;
      
      S.DATA.use.keyboard = 0;
      S.DATA.use.mouse = 0;
      S.DATA.use.eye = 0;
      
      % % % open the devices % % %
      Experiment.devices_open(S,DEVICES);
      
    end
    
    % % % Open devices and associated data % % %
    function data_open(S,dataType)
      
      % % % initialize the data storage structure % % %
      S.DATA.(dataType).table = table.empty();
      
      % % % the table for writing % % %
      S.DATA.(dataType).tableFull = table.empty();
      
      % % % information needed to build table % % %
      S.DATA.(dataType).columns = table.empty();
      
    end
    
    %% Run
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function run_subject(S,D)
      % record meta information
      S.META.didRun = 1;
      S.META.runTime = datestr(now);
      S.META.computername = getenv('computername');
      
      % run the block sequence script
      S.data_init(D.DEVICES);
      D.run(S);
      
      % record meta information
      S.META.didFinish = 1;
      S.META.finishTime = datestr(now);
      
    end
    
    %% SAVE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function save_subject(S)
      % save into the specified folder
      saveName = [pwd filesep 'DATA' filesep S.DATA.meta.fileHeader{:}];
      save(saveName,'S');
    end
  end
  
  methods (Static)
  end
end