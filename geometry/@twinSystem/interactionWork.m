function IW = interactionWork(tS,sigma, varargin)
%INTERACTION_WORK Summary of this function goes here
%   Detailed explanation goes here

% Interaction work with respect to a tension direction
% https://www.sciencedirect.com/science/article/pii/S0749641922002467?ref=pdf_download&fr=RR-2&rr=8c8ba49e3a256aa1#fig0002
% A positive value of IW means that the work is given by the external loading and the work could be consumed by the twin formation (Larcher et al., 2019). In contrast, a negative value of IW would mean that the material has transmitted energy to the machine, which is physically unrealistic.

if nargin == 1 || (isnumeric(sigma) && isempty(sigma))

  IW = S2FunHarmonic.quadrature(@(v) tS.interactionWork(v,varargin{:}),'bandwidth',4);

elseif isa(sigma,'stressTensor')
  sig = sigma.normalize;
  warning('interactionWork:stressTensor','Stress tensor input is normalized.')
  IW = sig : transpose(tS.displacementGradient);
  %trace(sig' .* tS.deformationTensor);
  % SF = dot_outer(r,eta1,'noSymmetry') .* dot_outer(r,k1,'noSymmetry');

elseif isa(sigma,'vector3d')

  r = sigma.normalize;
  sig = stressTensor.uniaxial(r);

  IW = sig : transpose(tS.deformationTensor); %trace(sig' .* tS.deformationTensor);
  % SF = dot_outer(r,eta1,'noSymmetry') .* dot_outer(r,k1,'noSymmetry');
end
