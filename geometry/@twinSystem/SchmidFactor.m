function SF = SchmidFactor(tS,sigma,varargin)
% compute the Schmid factor 
%
% Syntax
%
%   SF = SchmidFactor(sS,v)
%   SF = SchmidFactor(sS,sigma)
%   SF = SchmidFactor(sS,sigma,'relative')
%
% Input
%  tS - list of @twinSystem
%  v  - @vector3d list of tension direction
%  sigma - @stressTensor
%
% Output
%  SFfun - size(sS) x 1 list of @S2FunHarmonic
%  SF - size(sS) x size(sigma) matrix of Schmid factors
%   

eta1 = tS.eta1.normalize; %#ok<*PROPLC>
k1 = tS.k1.normalize;

% compute the relative Schmid factor by dividing by the critical resolved
% shear stress for every slip system
if check_option(varargin,'relative')
   eta1 = eta1./ tS.CRSS;
end

% % Schmid factor with respect to a tension direction
% if nargin == 1 || (isnumeric(sigma) && isempty(sigma))
% 
%   SF = S2FunHarmonic.quadrature(@(v) tS.SchmidFactor(v,varargin{:}),'bandwidth',4);
% probably something wrong in the subsref

if isa(sigma,'vector3d')
  
  r = sigma.normalize;
  SF = dot_outer(r,eta1,'noSymmetry') .* dot_outer(r,k1,'noSymmetry');
  
% Schmid factor with respect to a stress tensor
elseif isa(sigma,'stressTensor')
  
  % normalize the stress tensor
  % such that the resulting Schmid factor is always between [0, 0.5]
  EV = eig(sigma);
  sigma = sigma ./ reshape(EV(3,:)-EV(1,:),size(sigma));
  
  if isscalar(sigma)
    SF = EinsteinSum(sigma,[-1,-2],k1,-1,eta1,-2);
    SF = reshape(SF,size(tS));
  else
    SF = zeros(length(sigma),length(eta1));
  
    for i = 1:length(tS.eta1)
      SF(:,i) = EinsteinSum(sigma,[-1,-2],k1(i),-1,eta1(i),-2);
    end
  end
    
else
  
  error('Second argument should be either vector3d or stressTensor.')
  
end
end
