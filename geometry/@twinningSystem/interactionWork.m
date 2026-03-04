function IW = interactionWork(tS,sigma, varargin)
%INTERACTION_WORK Summary of this function goes here
%   Detailed explanation goes here

% Interaction work with respect to a tension direction
% https://www.sciencedirect.com/science/article/pii/S0749641922002467?ref=pdf_download&fr=RR-2&rr=8c8ba49e3a256aa1#fig0002
if nargin == 1 || (isnumeric(sigma) && isempty(sigma))
  
  IW = S2FunHarmonic.quadrature(@(v) tS.interactionWork(v,varargin{:}),'bandwidth',4);
    
elseif isa(sigma,'vector3d')
  
  r = sigma.normalize;
  sig = stressTensor.uniaxial(r);
  
  IW = sig : transpose(tS.deformationTensor); %trace(sig' .* tS.deformationTensor);
  % SF = dot_outer(r,eta1,'noSymmetry') .* dot_outer(r,k1,'noSymmetry');
end

