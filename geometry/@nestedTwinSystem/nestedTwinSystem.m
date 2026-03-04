classdef nestedTwinSystem
    properties
        primaryTwin    % Can be a twinSystem OR another nestedTwinSystem
        secondaryTwin  % twinSystem object
        CS             % Crystal symmetry (inherited from primary)
    end
    
    properties (Dependent)
        % These act like normal properties but are calculated on the fly,
        % allowing this class to mimic a standard twinSystem.
        k1             % Effective k1 in the parent reference frame
        eta1           % Effective eta1 in the parent reference frame
    end
    
    methods
        function nTS = nestedTwinSystem(tS_primary, tS_secondary)
            % Constructor
            nTS.primaryTwin = tS_primary;
            nTS.secondaryTwin = tS_secondary;
            nTS.CS = tS_primary.CS; 
        end
        
        function val = get.k1(nTS)
            % Rotates the secondary twin's habit plane back into the parent
            % coordinate system recursively.
            m1 = nTS.primaryTwin.parentTwinMisorientation();
            val = m1 * nTS.secondaryTwin.k1;
        end
        
        function val = get.eta1(nTS)
            % Rotates the secondary twin's shear direction back into the parent
            % coordinate system recursively.
            m1 = nTS.primaryTwin.parentTwinMisorientation();
            val = m1 * nTS.secondaryTwin.eta1;
        end
        
        function misori = parentTwinMisorientation(nTS)
            % Renamed to match twinSystem.m exactly.
            % Calculates the combined misorientation recursively.
            m1 = nTS.primaryTwin.parentTwinMisorientation();
            m2 = nTS.secondaryTwin.parentTwinMisorientation();
            
            % Creates a matrix of misorientations (Primary Variants x Secondary Variants)
            misori = m1(:) * m2(:).'; 
        end
        
        function defTensor = displacementGradient(nTS)
            % Renamed to match twinSystem.m exactly.
            % F_total = F_primary * F_secondary
            % H_total = F_total - Eye
            
            F1 = tensor.eye(nTS.primaryTwin.CS) + nTS.primaryTwin.displacementGradient();
            F2 = tensor.eye(nTS.secondaryTwin.CS) + nTS.secondaryTwin.displacementGradient();
            
            F_tot = F1 * F2; 
            defTensor = F_tot - tensor.eye(nTS.CS);
        end
        
        function SF = schmidFactor(nTS, stressTensor)
            % Calculates the Schmid factor using the effective parent elements.
            % Because k1 and eta1 are dependent properties, calling nTS.k1 
            % automatically triggers the rotation math.
            
            k1_vec = nTS.k1.xyz;
            eta1_vec = nTS.eta1.xyz;
            
            % SF = (eta1 * Stress * k1)
            SF = sum((eta1_vec * stressTensor) .* k1_vec, 2); 
        end

        function nTS_new = mtimes(obj1, obj2)
            % Overloads the * operator to create deeper nested twins.
            % This allows for syntax like: tertiaryTwin = doubleTwin * twin3;
            nTS_new = nestedTwinSystem(obj1, obj2);
        end
    end
end