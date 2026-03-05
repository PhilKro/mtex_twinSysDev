function tS = cat(dim, varargin)
% implement cat for twinSystem
%
% Syntax 
%   tS = cat(dim, tS1, tS2, tS3)
%
% Input
%  dim - dimension
%  varargin - @twinSystem objects
%
% Output
%  tS - @twinSystem

varargin(cellfun('isempty',varargin)) = [];
tS = varargin{1};

tSCS = cell(size(varargin)); tStwinType = tSCS; tSCRSS = tSCS; tSeta1= tSCS; tSeta2= tSCS; tSk1= tSCS; tSk2= tSCS; tSrotAxis= tSCS;
for i = 1:numel(varargin)
  tS2 = varargin{i};
  if ~isempty(tS2)
    tSeta1{i} = tS2.eta1;
    tSeta2{i} = tS2.eta2;
    tSk1{i} = tS2.k1;
    tSk2{i} = tS2.k2;
    tSrotAxis{i} = tS2.rotAxis;
    tSCRSS{i} = tS2.CRSS;
  end
end


tS.eta1 = cat(dim,tSeta1{:});
tS.eta2 = cat(dim,tSeta2{:});
tS.k1 = cat(dim,tSk1{:});
tS.k2 = cat(dim,tSk2{:});
tS.rotAxis = cat(dim,tSrotAxis{:});
tS.CS = cat(dim,tSCS{:});
tS.twinType = cat(dim,tStwinType{:});
tS.CRSS = cat(dim,tSCRSS{:});