function tS = repmat(tS,varargin)
% replicate twin system

tS.k1 = repmat(tS.k1,varargin{:});
tS.eta1 = repmat(tS.eta1,varargin{:});
tS.k2 = repmat(tS.k2,varargin{:});
tS.eta2 = repmat(tS.eta2,varargin{:});
tS.rotAxis = repmat(tS.rotAxis,varargin{:});
tS.CRSS = repmat(tS.CRSS,varargin{:});
tS.twinType = repmat(tS.twinType,varargin{:});
tS.variantId = repmat(tS.variantId,varargin{:});

if ~isempty(tS.parent)
    tS.parent = repmat(tS.parent,varargin{:});
end
end