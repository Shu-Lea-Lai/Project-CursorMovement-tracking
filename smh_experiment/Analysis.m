classdef Analysis < handle
  properties
    % % % name % % %
    % for user reference
    name = '';
    
    % % % datatable % % %
    DT = table();
    
    % % % type % % %
    % what type of analysis, controls behavior
    type = '';
    
    % % % factors % % %
    % which variables to include
    factors = struct('names',{''},'n',0);
    
    % % % levels % % %
    % subsets inside factors
    levels = struct('names',{''},'n',0);
    
    % % % covariates % % %
    % continuous variables
    cov = struct();
    
    
    % % % marginal-means % % %
    marginal = struct();
    
  end
  events
    factorNamesChange
  end
  methods
    function A = Analysis(DT)
      A.DT = DT;
      addlistener(A, 'factorNamesChange', @A.meta_updateFactorLevels);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% DESCRIPTIVE %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute Marginal Means for each combination of factors
    function calc_marginalMeans(A,dvName)
      %
      A.factors.table = smh_allComb(A.levels.n+1)-1;
      for rr = 1:size(A.factors.table,1)
        tableRow = A.factors.table(rr,:);
        % intitialize row logic
        rowLogic = true(size(A.DT,1),1);
        for cc = 1:size(A.factors.table,2)
          % pull out the relevant labels
          % column label in data table
          factorName = A.factors.names{cc};
          % which index of the label to pull
          factorLevel = tableRow(cc);
          
          % if the value is 0, we marginalize over everything, have logical
          % true for all rows. if this isn't 0, then grab only those rows
          % where there is a match
          
          if factorLevel == 0
            % all true
            columnLogic = true;
          else
            % what's the actual value of the factor at this level
            factorValue = A.levels.names{cc}(factorLevel);
            % what's the class
            factorClass = class(factorValue);
            switch factorClass
              case 'double'
                columnLogic = A.DT.(factorName) == factorValue;
              case 'cell'
                columnLogic = strcmp(A.DT.(factorName),factorValue);
            end
          end
          rowLogic = rowLogic & columnLogic;
        end
        % pull out the values from each row where the values were met
        vals = A.DT.(dvName)(rowLogic);
        A.marginal.n(rr) = length(vals);
        A.marginal.means(rr) = mean(vals);
        A.marginal.std(rr) = std(vals);
        A.marginal.se(rr) = A.marginal.std(rr) ./ sqrt(A.marginal.n(rr));
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% LISTENERS %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function set.factors(A,factors)
      A.factors = factors;
      A.factors.n = length(A.factors.names);
      A.meta_updateFactorLevels(factors)
    end
    
    function meta_updateFactorLevels(A,names)
      %       A.factors.names = names;
      A.levels = struct();
      A.levels.n = nan(1,A.factors.n);
      for ff = 1:A.factors.n
        factorName = A.factors.names{ff};
        A.levels.names(ff) = {unique(A.DT.(factorName))};
        A.levels.n(ff) = length(A.levels.names{ff});
      end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% MAIN EFFECTS %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Main effects functions
    
  end
  methods (Static)
  end
end