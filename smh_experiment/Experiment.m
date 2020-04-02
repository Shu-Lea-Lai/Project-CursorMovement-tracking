classdef Experiment < handle
  properties
    % % % FILE % % %
    FILE = struct();
    
    % % % META % % %
    META = struct();
    
    % % % SUBJECTS % % %
    SUBJECTS = table.empty();
    
    % % % DEVICES % % %
    DEVICES = struct();
    
  end
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% FILE MGMT %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% New
    function E = Experiment(experimentName)
      % construct an experiment object, which contains information about
      % the conditions, the subjects, the data, the analyses, and plotting
      
      % % % meta information % % %
      E.FILE.experimentName = experimentName;
      
      %       % % % add this folder to the path for now % % %
      addpath(experimentName);
      
      % see if this experiment already exists
      if E.experimentExists
        
        msgStr = ['Experiment ''' E.FILE.experimentName ...
          ''' already exists. Loading existing experiment...'];
        disp(msgStr);
        
        
        % load the existing experiment
        E = Experiment.load(E.FILE.experimentName);
        
      else
        E.FILE.isSaved = 0;
        E.FILE.lastSaved = '';
        
        
      end
    end
    %% Save
    function save(E)
      % save the current experiment
      
      % make sure the experiment name isn't empty (signifying error
      % somewhere)
      if ~isempty(E.FILE.experimentName)
        if ~E.experimentExists
          % if we can't find the intended save location
          % make sure the directory exists
          mkdir(E.FILE.experimentName)
        end
        save([E.defaultSaveLoc],'E')
        E.FILE.isSaved = 1;
        E.FILE.lastSaved = datestr(now);
        
      else
        errordlg('Experiment Name is empty, cannot save')
      end
    end
    %% Misc
    function foundExp = experimentExists(E)
      % looks to see if a folder with a name corresponding to
      % E.FILE.experimentName exists in the current working directory
      foundExp = exist(E.defaultSaveLoc,'file') > 0;
      
    end
    
    function fName = defaultSaveLoc(E)
      % default location to store an experiment
      fName = [E.FILE.experimentName filesep E.FILE.experimentName '_EXPERIMENT.mat'];
    end
    
    
    
    %% Classes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% RUNNING %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % create the initial class % % %
    % provides access to order information
    function C = class_create(E)
      C = [];
      eval(['C = ' E.FILE.experimentName ';']);
    end
    
    % % % initialize the design % % %
    % essentially generates all the information needed to actually run the
    % experiment such as creating the trials, etc...
    function C = class_init(E,S)
      C = E.class_create;
      C.init_Design(S);
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% SUBJECTS %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % load the datatable % % %
    function T = subs_tableLoad(E)
      dataName = [E.FILE.experimentName 'dataTable.txt'];
      % open it if it exists; otherwise create
      if exist(dataName)
        T = readtable(dataName,'filetype','text','DatetimeType','text');
      else
        T = array2table(cell(0,1),'VariableNames',{'DateTime'});
        T.SubjectNumber = nan(0,1);
        T.Order = nan(0,1);
        T.Booth = nan(0,1);
        T.Experimenter = cell(0,1);
      end
    end
    
    % % % save the datatable % % %
    function subs_tableSave(E,T)
      dataName = [E.FILE.experimentName 'dataTable.txt'];
      writetable(T,dataName);
    end
    
    
    % % % add a new subject % % %
    function S = subs_add(E)
      % load the existing table
      T = E.subs_tableLoad;
      % get the last subject in the table
      nSubs = size(T.SubjectNumber,1);
      boothNum = {'nan'};
      expName = {'EXPERIMENTER NAME'};
      if nSubs > 0
        % % % last sub info % % %
        lastSub = T(end,:);
        subNum = {num2str(lastSub.SubjectNumber + 1)};
      else
        subNum = {'101'};
      end
      
      orderNum = {num2str(E.subs_counterbalance)};
      % % % get inputs % % %
      
      % subject number
      in_subjectNumber = inputdlg('Subject Number:','Please enter the following',[1 50],subNum);
      in_orderNumber = inputdlg('Order:','Please enter the following',[1 50],orderNum);
      in_boothNumber = inputdlg('Booth Number:','Please enter the following',[1 50],boothNum);
      in_experimenterName = inputdlg('Experimenter Initials:','Please enter the following',[1 50],expName);
      
      % store and save into table
      T.DateTime{nSubs+1} = datestr(now);
      T.Experimenter(nSubs + 1) = in_experimenterName;
      T.SubjectNumber(nSubs+1) = str2num(in_subjectNumber{:});
      T.Order(nSubs + 1) = str2num(in_orderNumber{:});
      T.Booth(nSubs + 1) = str2num(in_boothNumber{:});
      
      
      % % % create a subject % % %
      S = Subject(E.FILE.experimentName);
      
      S.subjectNumber = T.SubjectNumber(end);
      S.config.order = T.Order(end);
      R = E.class_create;
      S.config.seed = R.info.orderSeeds(S.config.order);
      
      % if the experimenter inputs the subject number as 999, don't save
      % this person into the file
      if str2num(in_subjectNumber{:}) ~= 999
        E.subs_tableSave(T)
      end
    end
    
    % % % counterbalance % % %
    function [order,seed] = subs_counterbalance(E)
      R = E.class_create;
      T = E.subs_tableLoad;
      nSubs = size(T,1);
      % % % loop through % % %
      if nSubs > 0
        % how many subjects do we have in each order already?
        countOrders = zeros(1,R.info.nOrders);
        for ss = 1:length(T.SubjectNumber)
          thisOrder = T.Order(ss);
          countOrders(thisOrder) = countOrders(thisOrder) + 1;
        end
        % which has fewer?
        [~,least] = min(countOrders);
        % assign this order to the next subject
        order = least;
        seed = R.info.orderSeeds(order);
      else
        order = 1;
        seed = R.info.orderSeeds(1);
      end
      
    end
    
    %% Run Subject
    function subs_run(E,S)
      % seed the design
      D = E.class_init(S);
      S.run_subject(D);
      % thank and close
      Experiment.instru_thankyou(D.DEVICES);
      Experiment.devices_close
    end
    
    
  end
  
  methods(Static)
    %% FILES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% FILE MGMT %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load
    function E = load(name)
      E = load([name filesep name '_EXPERIMENT.mat']);
      E = E.E;
    end
    
    %% DEVICES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% DEVICES %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Open devices
    function devices_open(S,DEVICES)
      fn = fieldnames(DEVICES);
      nDevices = length(fn);
      
      % turn off the mouse by default
      HideCursor;
      
      for dd = 1:nDevices
        D = DEVICES.(fn{dd});
        [deviceName,isOutput] = Experiment.device_type(D);
        if isOutput
          S.data_open(deviceName);
          S.DATA.use.(deviceName) = 1;
        end
        
        % turn on the mouse if it's one of the devices
        if strcmp(deviceName,'mouse')
          ShowCursor(10);
        end
      end
    end
    
    function devices_close()
      % % % close the window % % %
      Device_Output_Display.close_all;
      % % % stop collecting kb % % %
      KbQueueStop
      % % % sound % % %
      try
        Device_Output_Sound.close_all
      end
    end
    
    % % % SCREEN % % %
    function info = init_device_window(opMode)
      switch opMode
        case 'full'
          info = Device_Output_Display.open_full;
        case {'demo','test'}
          info = Device_Output_Display.open_demo;
      end
    end
    % % % SOUND % % %
    function info = init_device_sound()
      info = [];
    end
    
    % % % DEVICE TYPE % % %
    function [deviceName,isOutput] = device_type(D)
      switch class(D)
        case 'Device_Input_Keyboard'
          deviceName = 'keyboard';
          isOutput = 1;
        case 'Device_Input_Mouse'
          deviceName = 'mouse';
          isOutput = 1;
        case 'Device_Input_Eyetracker'
          deviceName = 'eye';
          isOutput = 1;
        case 'Device_Output_Display'
          deviceName = 'display';
          isOutput = 0;
        case 'Device_Output_Sound'
          deviceName = 'sound';
          isOutput = 0;
      end
    end
    
    
    %% MISC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% MISC %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function globalKeyEsc = global_escape
      globalKeyEsc = KbName('esc');
      global globalOkay
      globalOkay = 1;
    end
    
    %% DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store the data
    function DATA = store_data(TABLE,S,mode)
      % % % add initial headers if missing % % %
      % if there aren't any headers, add them to the file
      tableNames = TABLE.Properties.VariableNames;
      DATA = S.DATA.(mode);
      
      if isempty(DATA.columns)
        DATA.columns = table(tableNames',[1:length(tableNames)]','VariableNames',{'Labels','Index'});
      end
      
      % % % update the table with new columns if necessary % % %
      % initialize the output
      printOut = table();
      for nn = 1:length(tableNames)
        cellMatch = strfind(DATA.columns.Labels,tableNames{nn});
        % find where it matches; must include the text and be the same
        % length
        logicalMatch = cellfun(@sum, cellMatch) & ...
          (cellfun(@length,DATA.columns.Labels) == length(tableNames{nn}));
        matchNum = find(logicalMatch);
        % if there wasn't a match found, add a column
        if isempty(matchNum)
          % add it to the info
          DATA.columns(end+1,:) = table(tableNames(nn),size(DATA.columns,1)+1);
          % then get the required index
          cellMatch = strfind(DATA.columns.Labels,tableNames{nn});
          logicalMatch = cellfun(@sum, cellMatch);
          matchNum = find(logicalMatch);
          % fill in the full table with nans here
          DATA.tableFull.(tableNames{nn}) = nan(size(DATA.tableFull,1),1);
        end
        % % % write the data % % %
        tableData = TABLE.(tableNames{nn});
        % store into the print out
        printOut.(tableNames{nn}) = tableData;
      end
      
      
      % % % update structure % % %
      printOutSafe = table_checkClasses(DATA.tableFull,printOut);
      printOutOrganized = table_arrangeColumns(DATA.columns,printOutSafe);
      DATA.tableFull(end+1,:) = printOutOrganized;
      
      
      % % % put back into subject % % %
      S.DATA.(mode) = DATA;
      
      function P = table_checkClasses(D,P)
        variableNames = D.Properties.VariableNames;
        nF = length(variableNames);
        for ff = 1:nF
          c1 = class(D.(variableNames{ff}));
          try
            c2 = class(P.(variableNames{ff}));
          catch
            switch c1
              case 'double'
                P.(variableNames{ff}) = nan;
              case 'cell'
                P.(variableNames{ff}) = {'nan'};
              case 'char'
                P.(variableNames{ff}) = 'nan';
            end
            c2 = class(P.(variableNames{ff}));
          end
          if ~strcmp(c1,c2)
            P.(variableNames{ff}) = cast(P.(variableNames{ff}),c1);
          end
        end
      end
      
      function T2 = table_arrangeColumns(C,T)
        T2 = table();
        nC = size(C,1);
        for cc = 1:nC
          cI = C.Index(cc);
          T2.(C.Labels{cI}) = T.(C.Labels{cI});
        end
        
      end
    end
    
    % Write the Data
    function DATA = write_data(S,mode)
      
      DATA = S.DATA.(mode);
      FID = DATA.fileID;
      COL = DATA.columns;
      % % % add headers % % %
      Experiment.write_headers(FID,COL);
      % format for printing out
      fmt = [repmat('%s ', 1, length(COL.Labels)) '\n'];
      % convert to a cell
      DATA.cell = table2cell(DATA.tableFull);
      
      for rr = 1:size(DATA.cell,1)
        tableRow = DATA.cell(rr,:);
        fprintf(FID,fmt,tableRow{:});
      end
    end
    
    
    % Write Column Headers
    function write_headers(FID,COL)
      % set the pointer to the beginning
      fseek(FID,0,-1);
      nC = size(COL,1);
      fprintf(FID, [repmat('%s ',1,nC) '\n'], COL.Labels{:});
      % move the pointer back to the end
      fseek(FID,0,1);
    end
    
    
    % % % concatenate data tables % % %
    % makes one large table that contains all of the sub-tables for a given
    % subject
    function data_combine(S)
      fn = fieldnames(S.DATA);
      tableFull = [];
      for ff = 1:length(fn)
        % don't try to combine 'meta' or 'use'
        if ~strcmp(fn{ff},'meta') && ~strcmp(fn{ff},'use') && isfield(S.DATA.(fn{ff}),'tableFull')
          tableFull = [tableFull,S.DATA.(fn{ff}).tableFull];
        end
      end
      % add a new field in DATA for this subject called 'all'
      S.DATA.all = tableFull;
    end
    
    %% Instructions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% INSTRUCTIONS %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % Intro Slide % % %
    function P = instru_welcome(S,name)
      welcomeText = Stimulus('text',...
        ['Welcome to the ' name ' Experiment'],... % str
        [0 0 0],... % color
        36,... % fontsize
        0.5,0.5); % posX, posY
      
      subText = Stimulus('text',...
        ['SUBJECT: ' num2str(S.subjectNumber)],...
        [0 0 0],... % color
        24,... % fontsize
        0.5,0.85); % posX, posY
      
      orderText = Stimulus('text',...
        ['Order: ' num2str(S.config.order)],...
        [0 0 0],... % color
        24,... % fontsize
        0.5,0.9); % posX, posY
      
      device = {'keyboard','button',KbName({'space','left','right'})};
      P = Phase([welcomeText,subText,orderText],nan,device);
    end
    
    % % % Thank You Slide % % %
    function instru_thankyou(D)
      thankText = Stimulus('text',...
        ['Thank you for participating!'],... % str
        [0 0 0],... % color
        36,... % fontsize
        0.5,0.25); % posX, posY
      
      experiText = Stimulus('text',...
        ['Please notify the experimenter that you are finished'],... % str
        [0 0 0],... % color
        36,... % fontsize
        0.5,0.5); % posX, posY
      
      device = {'keyboard','button',KbName({'space','left','right'})};
      P = Phase([thankText,experiText],60,device);
      
      P.run_Phase(D)
      
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function DT = analyze_import(varargin)
      % where are the files located
      dataFolder = [pwd filesep 'DATA' filesep];
      dataContents = dir([dataFolder '*.mat']);
      
      
      DT = table();
      % loop through and grab them all
      for dd = 1:length(dataContents)
        subjectData = load([dataFolder dataContents(dd).name]);
        S = subjectData.S;
        % condense the data
        Experiment.data_combine(S);
        % do some checks to make sure it's usable data
        if nargin > 0
          if size(S.DATA.all,1) == varargin{1}
            DT = [DT;S.DATA.all];
          end
        else
          DT = [DT;S.DATA.all];
        end
      end
    end
    
  end
end
