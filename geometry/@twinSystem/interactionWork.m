function IW = interactionWork(tS, sigma, varargin)
% compute the interaction work of a twin system under stress
%
% Syntax
%
%   IW = interactionWork(tS, v)
%   IW = interactionWork(tS, v, parentOri)
%   IW = interactionWork(tS, sigma)
%   IW = interactionWork(tS, sigma, parentOri)
%
% Input
%  tS - list of @twinSystem
%  v  - @vector3d (tension direction in specimen coords) or @Miller (in crystal coords)
%  sigma - @stressTensor (in specimen or crystal coords)
%  parentOri - @orientation (required if sigma or v is in specimen coordinates)
%
% Output
%  IW - interaction work

% Interaction work with respect to a tension direction
% https://www.sciencedirect.com/science/article/pii/S0749641922002467?ref=pdf_download&fr=RR-2&rr=8c8ba49e3a256aa1#fig0002
% A positive value of IW means that the work is given by the external loading and the work could be consumed by the twin formation (Larcher et al., 2019). In contrast, a negative value of IW would mean that the material has transmitted energy to the machine, which is physically unrealistic.

if nargin == 1 || (isnumeric(sigma) && isempty(sigma))
  IW = S2FunHarmonic.quadrature(@(v) tS.interactionWork(v,varargin{:}),'bandwidth',4);
  return;
end

% Check for optional parentOrientation in varargin
parentOri = orientation.empty;
for i = 1:length(varargin)
    if isa(varargin{i}, 'orientation')
        parentOri = varargin{i};
        break;
    end
end

% If a vector3d (or Miller) is supplied, instantiate a uniaxial stress tensor
if isa(sigma, 'vector3d')
    if isa(sigma, 'Miller')
        % It is in crystal coordinates, check if CS matches
        if ~strcmp(sigma.CS.mineral, tS.CS.mineral)
            warning('interactionWork:CSMismatch', 'The crystal symmetry of the Miller direction does not match the twin system.');
        end
        % Instantiate uniaxial stress state along this crystal direction
        sigma = stressTensor.uniaxial(sigma);
    else
        % It is a standard vector3d in specimen coordinates
        if isempty(parentOri)
            error('interactionWork:MissingOrientation', 'A parent orientation must be supplied when providing a tension direction in specimen coordinates.');
        end
        % Transform the uniaxial stress tensor to crystal coordinates
        sigma = inv(parentOri) * stressTensor.uniaxial(sigma);
    end
end

% At this point, sigma is guaranteed to be a stressTensor
if isa(sigma, 'stressTensor')
    
    % Check if the stress tensor is in specimen or crystal coordinates
    if isa(sigma.CS, 'specimenSymmetry')
        if isempty(parentOri)
            error('interactionWork:MissingOrientation', 'A parent orientation must be supplied when providing a stress tensor in specimen coordinates.');
        end
        % Transform to crystal coordinates
        sigma = inv(parentOri) * sigma;
    elseif isa(sigma.CS, 'crystalSymmetry')
        if ~strcmp(sigma.CS.mineral, tS.CS.mineral)
            warning('interactionWork:CSMismatch', 'The crystal symmetry of the stress tensor does not match the twin system.');
        end
    end
    
    % sig = sigma.normalize;
    % warning('interactionWork:stressTensor','Stress tensor input is normalized.')
    sig = sigma;
    IW = sig : transpose(tS.displacementGradient); 
    
else
    error('interactionWork:InvalidInput', 'Second argument should be either vector3d, Miller, or stressTensor.')
end
end
