function varargout = size(tS,varargin)
% overloads size

[varargout{1:nargout}] = size(tS.eta1.x,varargin{:});