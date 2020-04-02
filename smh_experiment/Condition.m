classdef Condition < handle
  properties
    settings = struct();
  end
  methods
    function C = Condition(properties)
      fn = fieldnames(properties);
      for ff = 1:length(fn)
        C.settings.(fn{ff}) = properties.(fn{ff});
      end
    end
  end
  methods (Static)
  end
end