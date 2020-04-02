% function [currentState, vStore, logP] = hmmviterbi_multi2(seq,tr,e,varargin)

[seq, ~] = hmmgenerate(100,tr,e);
tr = [0.95,0.05;
          0.10,0.90];

e = [1/6,  1/6,  1/6,  1/6,  1/6,  1/6;
         1/10, 1/10, 1/10, 1/10, 1/10, 1/2;];
     
     
numStates = size(tr,1);
checkTr = size(tr,2);
if checkTr ~= numStates
  error(message('stats:hmmviterbi:BadTransitions'));
end

% number of rows of e must be same as number of states

checkE = size(e,1);
if checkE ~= numStates
  error(message('stats:hmmviterbi:InputSizeMismatch'));
end

numEmissions = size(e,2);
customStatenames = false;

% deal with options
if nargin > 3
  okargs = {'symbols','statenames'};
  [symbols,statenames] = ...
    internal.stats.parseArgs(okargs, {[] []}, varargin{:});
  
  if ~isempty(symbols)
    numSymbolNames = numel(symbols);
    if ~isvector(symbols) || numSymbolNames ~= numEmissions
      error(message('stats:hmmviterbi:BadSymbols'));
    end
    [~, seq]  = ismember(seq,symbols);
    if any(seq(:)==0)
      error(message('stats:hmmviterbi:MissingSymbol'));
    end
  end
  if ~isempty(statenames)
    numStateNames = length(statenames);
    if numStateNames ~= numStates
      error(message('stats:hmmviterbi:BadStateNames'));
    end
    customStatenames = true;
  end
end


% work in log space to avoid numerical issues
L = length(seq);
if any(seq(:)<1) || any(seq(:)~=round(seq(:))) || any(seq(:)>numEmissions)
  error(message('stats:hmmviterbi:BadSequence', numEmissions));
end
currentState = zeros(1,L);
if L == 0
  return
end
logTR = log(tr);
logE = log(e);

% allocate space
pTR = zeros(numStates,L);
% assumption is that model is in state 1 at step 0
v = -Inf(numStates,1);
v(1,1) = 0;
vOld = v;

vStore = nan(size(pTR));

% loop through the model
for count = 1:L
  for state = 1:numStates
    % for each state we calculate
    % v(state) = e(state,seq(count))* max_k(vOld(:)*tr(k,state));
    bestVal = -inf;
    bestPTR = 0;
    % use a loop to avoid lots of calls to max
    for inner = 1:numStates
      val = vOld(inner) + logTR(inner,state);
      if val > bestVal
        bestVal = val;
        bestPTR = inner;
      end
    end
    % save the best transition information for later backtracking
    pTR(state,count) = bestPTR;
    %%% update v - changed from default
    vCalc = 0;
    for meas = 1:size(logE,3)
      vCalc = vCalc + logE(state,seq(state,count,meas),meas);      
    end
    v(state) = vCalc + bestVal;
  end
  vStore(:,count) = v-vOld; 
  vOld = v;  
end

% decide which of the final states is post probable
[logP, finalState] = max(v);

% Now back trace through the model
currentState(L) = finalState;
for count = L-1:-1:1
  currentState(count) = pTR(currentState(count+1),count+1);
  if currentState(count) == 0
    error(message('stats:hmmviterbi:ZeroTransitionProbability', currentState( count + 1 )));
  end
end
if customStatenames
  currentState = reshape(statenames(currentState),1,L);
end


