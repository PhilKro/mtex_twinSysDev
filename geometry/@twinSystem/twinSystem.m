classdef twinSystem
    properties
        CRSS     % Critical resolved shear stress
        eta1     % Shear direction
        eta2     % Conjugate shear direction
        k1       % Twin plane
        k2       % Conjugate twin plane
        rotAxis  % Zone axis
        twinType % Type of twin: 0:'Compound', 1:'Type I' or 2:'Type II'
        parent   % Parent object (orientation, grain2d, or twinSystem)
        variantId % Variant ID
        
    end

    properties (Dependent = true)
        CS       % Crystal symmetry
        isSymmetrised
    end

    methods
        function tS = twinSystem(rotAxis, k1, eta1, k2, eta2, CRSS, type, variantId)
            if nargin > 0
                tS.rotAxis = rotAxis;
                tS.k1 = k1;
                tS.eta1 = eta1;
                tS.k2 = k2;
                tS.eta2 = eta2;                

                % Set default twin type
                if nargin < 7 || isempty(type)
                    tS.twinType = 1;
                else
                    tS.twinType = type;
                end

                % Handle CRSS
                if nargin >= 6 && ~isempty(CRSS)
                    tS.CRSS = CRSS;
                else
                    tS.CRSS = 1; % Default value
                end

                % Handle variantId
                if nargin < 8 || isempty(variantId)
                    tS.variantId = ones(size(tS.k1));
                else
                    tS.variantId = variantId;
                end

                % Automatically set to compound if S (shear plane) is a mirror plane
                isC = tS.isCompound;
                if any(isC)
                    tS.twinType(isC) = 0;
                    %% TODO CHECK FOR ALLOWED PLANES AND DIRECTIONS FOR THE TWIN ELEMENTS BASED ON THE TWIN TYPE+
                end
            end
        end

        function CS = get.CS(tS)
        if isa(tS.eta1,'Miller')
            CS = tS.eta1.CS;
        else
            CS = specimenSymmetry.default;
        end
        end
        
        function out = get.isSymmetrised(tS)
        if length(tS)<2
            out = false;
        else
            out = eq(tS.subSet(1),tS.subSet(2));
        end
        end

        
        function [variantIndex] = variantDetermination2(tS, misOrientation, deviationInDegree)
            warning('twinSystem:variantDetermination2:deprecated', 'variantDetermination2 has to be deleted, also purged from subsref.')
            tSvariants = tS.symmetrise('antipodal');
            variantPTM = tSvariants.parentTwinMisorientation();
            if iscolumn(variantPTM)
                if iscolumn(misOrientation)
                    misOrientation = misOrientation.';
                end
            else
                if isrow(misOrientation)
                    misOrientation = misOrientation.'; 
                end
            end
            variantMatch_angle = angle(misOrientation, variantPTM, 'noSym1') / degree;
            variantMatch_logic = variantMatch_angle == min(variantMatch_angle, [], 2) & variantMatch_angle < deviationInDegree;
            
            variantIndex = variantMatch_logic * transpose(1:length(tSvariants));
        end
        

        function isC = isCompound(tS)
            % A twin is compound if the plane of shear S is a mirror plane.
            % The plane of shear S is the cross product of eta1 and k1.
            S = cross(tS.eta1, tS.k1);
            isC = S.isMirrorPlane;
        end

        function misori = parentTwinMisorientation(tS)
            warning('twinSystem:parentTwinMisorientation:deprecated', 'parentTwinMisorientation does not handle recursive twin objects yet.')
            % The misorientation depends on the twin type.
            misori = orientation.nan(tS.CS, tS.CS, size(tS));

            % Type I (1) and Compound (0)
            isType1OrCompound = tS.twinType == 1 | tS.twinType == 0;
            if any(isType1OrCompound)
                % Type I twin is defined by a reflection in the K1 plane.
                % For centrosymmetric crystals, this is equivalent to a 180-degree rotation about the normal to K1.
                warning('twinSystem:parentTwinMisorientation:type1Assumption', ...
                    'The current implementation works with an improper rotation here. This might not work for the orientation analysis. Maybe the proper rotation should be implemented here aswell/instead?');
                misori(isType1OrCompound) = -orientation.byAxisAngle(tS.k1(isType1OrCompound), pi);
            end
            % Type II (2)
            isType2 = tS.twinType == 2;
            if any(isType2)
                % Type II twin is defined by a 180-degree rotation about eta1.
                misori(isType2) = orientation.byAxisAngle(tS.eta1(isType2), pi);
            end
        end

        function shear = shearMagnitude(tS)
            warning('twinSystem:shearMagnitude:deprecated', 'shearMagnitude does not handle recursive twin objects yet. Might anyways be redundand with the shear method.')
            shear = sqrt(4 * ((dot(tS.eta2.normalize, tS.k1.normalize).^-2) - 1));
        end
        
        function [shearvector, shearMagnitude] = shear(tS)
            warning('twinSystem:shear:deprecated', 'Shearmagnitude does not handle recursive twin objects yet.')
            shearvector = 2 * (tS.k1.normalize - dot(tS.eta2.normalize, tS.k1.normalize).^-1 * tS.eta2.normalize);
            shearvector.dispStyle = 'UVTW';
            shearMagnitude = norm(shearvector);
        end
        
        function defTensor = displacementGradient(tS)
            warning('twinSystem:displacementGradient:deprecated', 'Hopefully displacementGradient handles recursive twin objects correctly.')
            % 1. Calculate the local deformation tensor in the parent crystal reference frame
            defTensor = tS.shearMagnitude .* dyad(tS.eta1, tS.k1).normalize;

            % 2. If no parent is attached, return the tensor in crystal coordinates
            if isempty(tS.parent) || all(cellfun('isempty', {tS.parent}))
                return;
            end

            % 3. Retrieve the root parent to check if it is tied to specimen coordinates
            superP = tS.superParent();

            if isa(superP, 'orientation') || isa(superP, 'grain2d')
                parents = [tS.parent];

                % Extract the specimen orientation of the immediate parent
                if isa(parents, 'twinSystem')
                    pOri = parents.orientation();
                elseif isa(parents, 'grain2d')
                    pOri = parents.meanOrientation;
                elseif isa(parents, 'orientation')
                    pOri = parents;
                else
                    return;
                end

                % Ensure the parent orientation array matches the dimensions of the twin array
                pOri = reshape(pOri, size(defTensor));

                % 4. Transform the deformation tensor from crystal to specimen coordinates
                defTensor = pOri * defTensor;
            end
        end
        
        function correspondanceMatrix = correspondanceMatrix(tS)
            warning('twinSystem:correspondanceMatrix:deprecated', 'correspondanceMatrix does not handle recursive twin objects yet.')
            correspondanceMatrix = -tensor.eye(tS.CS) + 2 * dyad(tS.eta2, tS.k1).normalize;
        end

        function nrVariants = numberOfVariants(tS)
            nrVariants = length(tS.symmetrise('antipodal'));
        end

        function display(tS, varargin) %#ok<DISPLAY>
            displayClass(tS, inputname(1), varargin{:}, 'moreInfo', char(tS.CS, 'compact'));

            if length(tS) > 24, disp([' CRSS: ' xnum2str(unique(tS.CRSS))]); end
            if length(tS) > 1, disp([' size: ' size2str(tS.CRSS)]); end % Replaced tS.b with tS.eta1 here

            disp(' ');

            if length(tS) <= 45 && ~isempty(tS)
                dispData(tS)
            elseif ~getMTEXpref('generatingHelpMode')
                disp(' ')
                setappdata(0, 'dispLastData', @() dispData(tS));
                disp('  <a href="matlab:feval(getappdata(0,''dispLastData''))">display all coordinates</a>')
                disp(' ')
            end
        end

        function n = char(tS)
            reta1 = char(round(tS.eta1));
            rk1 = char(round(tS.k1));
            if isscalar(tS)
                switch tS.twinType
                    case 0
                        typeName = "Compound";
                    case 1
                        typeName = "Type I";
                    case 2
                        typeName = "Type II";
                    otherwise
                        typeName = "Unknown";
                end
                n = strcat(strtrim(rk1), " ", strtrim(reta1), " (", typeName, ")");
            else
                n = strcat(rk1, " ", reta1);
            end
        end

        function n = filename_char(tS)
            reta1 = char(round(tS.eta1));
            rk1 = char(round(tS.k1));
            n = strcat("p", rk1(2:end-1), "d", reta1(2:end-1));
            n = strrep(n, '-', '_');
        end

        function dispData(tS)
            if isa(tS.CS, 'crystalSymmetry')
                typeNames = strings(length(tS.eta1),1);
                typeNames(tS.twinType == 0) = "Compound";
                typeNames(tS.twinType == 1) = "Type I";
                typeNames(tS.twinType == 2) = "Type II";

                if tS.eta1.lattice.isTriHex
                    reta1 = round(tS.eta1);
                    rk1 = round(tS.k1);
                    d = [reta1.UVTW rk1.hkil];
                    d(abs(d) < 1e-10) = 0;
                    dataCell = num2cell([d, reshape(tS.CRSS, [], 1), reshape(tS.variantId, [], 1)]);
                    fullRow = [dataCell, cellstr(typeNames)];
                    fullRow = fullRow.'; % Transpose for column-major printing
                    fprintf('eta_1 U   V   T   W | K_1 H   K   I   L   CRSS  Variant   Type \n');
                    fprintf('%7.0f %3.0f %3.0f %3.0f |  %3.0f %3.0f %3.0f %3.0f %6.2f %8d   %s\n', fullRow{:});
                    % cprintf([d, reshape(tS.CRSS, [], 1), char({tS.twinType})], '-L', '  ', '-Lc', {'eta_1 U' 'V' 'T' 'W' '| K_1 H' 'K' 'I' 'L' 'CRSS' 'Type'});
                else
                    d = [tS.eta1.uvw tS.k1.hkl];
                    d(abs(d) < 1e-10) = 0; 
                    numericData = [d, reshape(tS.CRSS, [], 1), reshape(tS.variantId, [], 1)];
                    dataCell = [num2cell(numericData), cellstr(typeNames)];
                    cprintf(dataCell, '-L', '  ', '-Lc', {'eta_1 u' 'v' 'w' '| K_1 h' 'k' 'l' 'CRSS' 'Variant' 'Type'});
                    % cprintf([d, reshape(tS.CRSS, [], 1), char({tS.twinType})], '-L', '  ', '-Lc', {'eta_1 u' 'v' 'w' '| K_1 h' 'k' 'l' 'CRSS' 'Type'});
                end
            else
                d = round(100 * [tS.eta1.xyz tS.k1.xyz]) ./ 100;
                d(abs(d) < 1e-10) = 0;
                cprintf(d, '-L', '  ', '-Lc', {'x' 'y' 'z' ' |   x' 'y' 'z'});
            end
        end

        function tS = mtimes(A, B)
            %MTIMES implement multiplication for twinSystem
            %
            % Syntax
            %   tS_nested = tS1 * tS2
            %   tS_active = ori * tS
            %   tS_active = grains * tS
            
            if isa(A, 'twinSystem') && isa(B, 'twinSystem')
                % Nesting: A is the parent of B
                parentObj = A(:);
                childObj = B(:);
                
                nParent = numel(parentObj.eta1);
                nChild = numel(childObj.eta1);
                
                % Expand child (inner loop) -> [C1; C2; C1; C2]
                tS = repmat(childObj, nParent, 1);
                
                % Expand parent (outer loop) -> [P1; P1; P2; P2]
                % Use indexing for robustness with object arrays
                idx = kron(1:nParent, ones(1, nChild)).';
                parentsExpanded = parentObj.subSet(idx);
                
                % Assign parents
                if nParent * nChild > 0
                    tS.parent = parentsExpanded;
                end
                
            elseif (isa(A, 'orientation') || isa(A, 'grain2d')) && isa(B, 'twinSystem')
                % Activation: A is the root parent
                tS = activate(B, A);
            else
                tS = builtin('mtimes', A, B);
                return;
            end
        end

        function tS = activate(tS, parentObj)
            % ACTIVATE Attach a parent (orientation/grain) to the twin system
            %
            % Syntax
            %   tS = activate(tS, parentObj)
           
            if ~ eq(parentObj.CS, tS.CS)
                error('twinSystem:activate:incompatibleCS', ...
                    'The crystal symmetry of the parent orientation/grain must match that of the twin system.');
            end

            if ~isa(parentObj, 'grain2d')
                parentObj = parentObj(:);
            end
            childObj = tS(:);
            
            nParent = numel(parentObj);
            if isa(parentObj, 'twinSystem')
                nParent = numel(parentObj.eta1);
            end

            nChild = numel(childObj);
            if isa(childObj, 'twinSystem')
                nChild = numel(childObj.eta1);
            end
            
            % Expand child (inner loop) -> [C1; C2; C1; C2]
            tS = repmat(childObj, nParent, 1);
            
            % Expand parent (outer loop) -> [P1; P1; P2; P2]
            idx = kron(1:nParent, ones(1, nChild)).';
            parentsExpanded = parentObj(idx);
            
            % Assign parents
            if nParent * nChild > 0
                tS.parent = parentsExpanded;
            end
        end

        function ori = orientation(tS)
            % ORIENTATION Calculate the twin orientation in specimen coordinates
            
            if isempty(tS),warning('twinSystem:orientation:empty', 'Your twin system is empty.'); ori = []; return; end
            
            % Check if parents are defined (handle object arrays safely)
            parentsCell = {tS.parent};
            if all(cellfun('isempty', parentsCell))
                warning('twinSystem:orientation:noParent', ...
                    'Parent orientation is not defined for any twin system in the array. Returning NaN orientations.');
                ori = orientation.nan(tS.CS,size(tS)); 
                return; 
            end

            % Recursively resolve parent orientation
            % This concatenation will fail if parents are of mixed types (e.g. grain and orientation)
            % which is desired behavior as we can't process mixed types easily.
            parents = [tS.parent];
            
            % PROBOABLY A PROBLEM WITH twinSystem as parent
            if isa(parents, 'twinSystem')
                numparents = numel(parents.eta1);
            else
                numparents = numel(parents);
            end

            if numparents ~= numel(tS.eta1)
                error('twinSystem:orientation:missingParent', ...
                    'Parent orientation is not defined for all twin systems in the array.');
            end

            if isa(parents, 'twinSystem')
                pOri = orientation(parents);
            elseif isa(parents, 'grain2d')
                pOri = parents.meanOrientation;
            else
                pOri = parents;
            end
            
            % Reshape pOri to match tS dimensions
            pOri = reshape(pOri, size(tS));

            % Apply misorientation (Parent -> Twin)
            ori = pOri .* tS.parentTwinMisorientation;
        end

        function [superP, n] = superParent(tS)
            % SUPERPARENT Find the highest level parent of the given object
            %
            % Syntax
            %   [superP, n] = superParent(tS)
            %
            % Output
            %   superP - The root parent object (orientation, grain2d, or twinSystem)
            %   n      - Number of layers traversed
            
            superP = tS;
            n = 0;
            
            while isa(superP, 'twinSystem') && ~isempty(superP)
                % Check if parents exist
                parentsCell = {superP.parent};
                
                % If any element is missing a parent, we stop to avoid returning
                % a mixed list or erroring on concatenation.
                if any(cellfun('isempty', parentsCell))
                    break;
                end
                
                % Move up one level
                try
                    nextParents = [superP.parent];
                catch
                    warning('twinSystem:superParent:mixedParents', ...
                        'Could not concatenate parents. They might be of mixed types.');
                    break;
                end
                
                superP = nextParents;
                n = n + 1;
            end
        end

        function n = length(tS)
            n = length(tS.k1);
        end

        % function n = numel(tS,varargin)
        %uncommenting this breaks the logic for some reason if you have a
        %nx1 twin system tS.parent calls this and then it fails in subsref
        %whatevs
        %     n = numel(tS.k1,varargin{:});
        % end

        function varargout = size(tS,varargin)
            [varargout{1:nargout}] = size(tS.k1,varargin{:});
        end

        function l = eq(tS1,tS2)
            if numel(tS1) > 1 || numel(tS2) > 1
                 if numel(tS1) == 1, tS1 = repmat(tS1,size(tS2)); end
                 if numel(tS2) == 1, tS2 = repmat(tS2,size(tS1)); end
                 l = false(size(tS1));
                 for i=1:numel(tS1)
                     l(i) = eq(tS1(i),tS2(i));
                 end
                 return
            end
            
            if tS1.CS ~= tS2.CS, l = false; return; end
            
            %% CHECK IF CORRECT
            warning('twinSystem:eq:correctness','TODO: check if the implementation of twinSystem/eq is correct.')
            ops = tS1.CS.quaternion;
            l = any( (ops * tS1.k1) == tS2.k1 & (ops * tS1.eta1) == tS2.eta1 );
        end

        function l = ne(tS1,tS2)
            l = ~eq(tS1,tS2);
        end
    end

    methods (Static = true)
        function tS = byK1eta2(k1, eta2, type)
            if nargin < 3, type = 1; end
            rotAxis = cross(k1, eta2);
            eta1 = cross(k1, rotAxis);
            k2 = cross(eta2, rotAxis);
            tS = twinSystem(rotAxis, k1, eta1, k2, eta2, 1, type);
            
            shearvec = 2 * (tS.k1.normalize - dot(tS.eta2.normalize, tS.k1.normalize).^-1 * tS.eta2.normalize);
            assert(dot(shearvec, tS.eta1, 'noSymmetry') > 0, "Something went wrong in the definition of the twin system: eta1 has the wrong sense. Sense of shear will not be correct")
        end

        function tS = byK2eta1(k2, eta1, type)
            if nargin < 3, type = 1; end
            rotAxis = cross(eta1, k2);
            rotAxis.dispStyle = 'UVTW';
            k1 = cross(rotAxis, eta1);
            k1.dispStyle = 'hkil';
            eta2 = cross(rotAxis, k2);
            eta2.dispStyle = 'UVTW';
            tS = twinSystem(rotAxis, k1, eta1, k2, eta2, 1, type);
            
            shearvec = 2 * (tS.k1.normalize - dot(tS.eta2.normalize, tS.k1.normalize).^-1 * tS.eta2.normalize);
            assert(dot(shearvec, tS.eta1, 'noSymmetry') > 0, "Something went wrong in the definition of the twin system: eta1 has the wrong sense. Sense of shear will not be correct")
        end

        function tS = byK1rotAxis(k1, rotAxis, type)
            if nargin < 3, type = 2; end
            if dot(k1, rotAxis, 'noSymmetry') < 1.00e-012
                eta1 = cross(k1, rotAxis);
                eta1.dispStyle = 'UVTW';
                k2 = eta1;
                k2.dispStyle = 'hkil';
                eta2 = k1;
                eta2.dispStyle = 'UVTW';
                disp(append(newline, 'Chose a set of values for K2 and eta2!'))
                tS = twinSystem(rotAxis, k1, eta1, k2, eta2, 1, type);
            else
                error("k1 and rotAxis have to be orthogonal")
            end
        end

        function tS = hexagonal_110K(K, CS)
            if CS.lattice.isTriHex
                tS = arrayfun(@(x_int) twinSystem.byK1rotAxis(Miller({1, 0, -1, x_int}, CS), Miller({1, -2, 1, 0}, CS, 'UVTW'), 2), K);
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_211K(K, CS)
            if CS.lattice.isTriHex
                tS = arrayfun(@(x_int) twinSystem.byK1rotAxis(Miller({1, 1, -2, x_int}, CS), Miller({1, -1, 0, 0}, CS, 'UVTW'), 2), K);
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1012_1012(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.byK2eta1(Miller(-1, 0, 1, 2, CS), Miller(-1, 0, 1, 1, CS, 'UVTW'), 1);
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1122_0001(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.byK2eta1(Miller(0, 0, 0, 1, CS), Miller(-1, -1, 2, 3, CS, 'UVTW'), 1);
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1101_1013(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.byK2eta1(Miller(1, 0, -1, -3, CS), Miller(1, 0, -1, -2, CS, 'UVTW'), 1);
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1122_1124(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.byK2eta1(Miller(1, 1, -2, -4, CS), Miller(1, 1, -2, -3, CS, 'UVTW'), 1);
            else
                disp("Provide a hexagonal CS for this method")
            end
        end
        
        function tS = hexagonal_1124_1122(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.byK2eta1(Miller(1, 1, -2, -2, CS), Miller(2, 2, -4, -3, CS, 'UVTW'), 1);
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1121_0001(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.byK2eta1(Miller(0, 0, 0, 1, CS), Miller(-1, -1, 2, 6, CS, 'UVTW'), 1);
            else
                disp("Provide a hexagonal CS for this method")
            end
        end
    end
end