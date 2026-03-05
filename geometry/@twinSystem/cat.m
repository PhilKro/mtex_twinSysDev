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

tSCRSS = cell(size(varargin)); tStwinType = tSCRSS; tSeta1= tSCRSS; tSeta2= tSCRSS; tSk1= tSCRSS; tSk2= tSCRSS; tSrotAxis= tSCRSS; tSparent = tSCRSS;
for i = 1:numel(varargin)
  tS2 = varargin{i};
  if ~isempty(tS2)
    tSeta1{i} = tS2.eta1;
    tSeta2{i} = tS2.eta2;
    tSk1{i} = tS2.k1;
    tSk2{i} = tS2.k2;
    tSrotAxis{i} = tS2.rotAxis;
    tStwinType{i} = tS2.twinType;
    tSCRSS{i} = tS2.CRSS;
    tSparent{i} = tS2.parent;
  end
end


tS.eta1 = cat(dim,tSeta1{:});
tS.eta2 = cat(dim,tSeta2{:});
tS.k1 = cat(dim,tSk1{:});
tS.k2 = cat(dim,tSk2{:});
tS.rotAxis = cat(dim,tSrotAxis{:});
tS.twinType = cat(dim,tStwinType{:});
tS.CRSS = cat(dim,tSCRSS{:});
tS.parent = cat(dim,tSparent{:});