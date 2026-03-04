function tS = cat(dim,varargin)
% implement cat for slipSystem
%
% Syntax 
%   tS = cat(dim,tS1,tS2,tS3)
%
% Input
%  dim - dimension
%  tS1, tS2, tS3 - @slipSystem
%
% Output
%  tS - @twinningSystem
%
% See also
%TODO
% twinningSystem/horzcat, twinningSystem/vertcat

% remove emtpy arguments
varargin(cellfun('isempty',varargin)) = [];
tS = varargin{1};

tSb = cell(size(varargin)); tSn = tSb; tSCRSS = tSb; tSeta1= tSb; tSeta2= tSb; tSk1= tSb; tSk2= tSb; tSrotAxis= tSb;
for i = 1:numel(varargin)
  tS2 = varargin{i};
  if ~isempty(tS2)
    tSb{i} = tS2.b;
    tSn{i} = tS2.n;
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
tS.b = cat(dim,tSb{:});
tS.n = cat(dim,tSn{:});
tS.CRSS = cat(dim,tSCRSS{:});