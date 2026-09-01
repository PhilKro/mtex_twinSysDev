function SF = SchmidFactor(tS, sigma, varargin)
% compute the Schmid factor 
%
% Syntax
%
%   SF = SchmidFactor(tS, v)
%   SF = SchmidFactor(tS, v, parentOri)
%   SF = SchmidFactor(tS, sigma)
%   SF = SchmidFactor(tS, sigma, parentOri)
%   SF = SchmidFactor(tS, sigma, 'relative')
%
% Input
%  tS - list of @twinSystem
%  v  - @vector3d (tension direction in specimen coords) or @Miller (in crystal coords)
%  sigma - @stressTensor (in specimen or crystal coords)
%  parentOri - @orientation (required if sigma or v is in specimen coordinates)
%
% Options
%  'relative' - compute the relative Schmid factor by dividing the Schmid factor by the critical resolved shear stress (CRSS) of each twin system
%
% Output
%  SF - size(tS) x size(sigma) matrix of Schmid factors
%   

eta1 = tS.eta1.normalize; %#ok<*PROPLC>
k1 = tS.k1.normalize;

% Check for optional parentOrientation in varargin
parentOri = orientation.empty;
for i = 1:length(varargin)
    if isa(varargin{i}, 'orientation')
        parentOri = varargin{i};
        break;
    end
end

% Compute the relative Schmid factor by dividing by the critical resolved
% shear stress for every slip system
if check_option(varargin,'relative')
   eta1 = eta1./ tS.CRSS;
end

% If a vector3d (or Miller) is supplied, instantiate a uniaxial stress tensor
if isa(sigma, 'vector3d')
    if isa(sigma, 'Miller')
        % It is in crystal coordinates, check if CS matches
        if ~strcmp(sigma.CS.mineral, tS.CS.mineral)
            warning('SchmidFactor:CSMismatch', 'The crystal symmetry of the Miller direction does not match the twin system.');
        end
        % Instantiate uniaxial stress state along this crystal direction
        sigma = stressTensor.uniaxial(sigma);
    else
        % It is a standard vector3d in specimen coordinates
        if isempty(parentOri)
            error('SchmidFactor:MissingOrientation', 'A parent orientation must be supplied when providing a tension direction in specimen coordinates.');
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
            error('SchmidFactor:MissingOrientation', 'A parent orientation must be supplied when providing a stress tensor in specimen coordinates.');
        end
        % Transform to crystal coordinates
        sigma = inv(parentOri) * sigma;
    elseif isa(sigma.CS, 'crystalSymmetry')
        if ~strcmp(sigma.CS.mineral, tS.CS.mineral)
            warning('SchmidFactor:CSMismatch', 'The crystal symmetry of the stress tensor does not match the twin system.');
        end
    end
    
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
    error('SchmidFactor:InvalidInput', 'Second argument should be either vector3d, Miller, or stressTensor.')
end
end
