classdef twinSystem
    properties
        CRSS     % Critical resolved shear stress
        eta1     % Shear direction
        eta2     % Conjugate shear direction
        k1       % Twin plane
        k2       % Conjugate twin plane
        rotAxis  % Zone axis
        twinType % Type of twin: 0:'Compound', 1:'Type I' or 2:'Type II'
        
    end

    properties (Dependent = true)
        CS       % Crystal symmetry
        isSymmetrised
    end

    methods
        function tS = twinSystem(rotAxis, k1, eta1, k2, eta2, CRSS, type)
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

        function [isTwinning, variantIndex, tSvariants] = variantDetermination(tS, misOrientation, deviationInDegree)
            tSvariants = tS.symmetrise('antipodal');
            variantMatch = angle(misOrientation, tSvariants.parentTwinMisorientation(), 'noSym1') / degree;
            variantMatch_logic = variantMatch == min(variantMatch, [], 2) & variantMatch < deviationInDegree;
            isTwinning = logical(variantMatch_logic * ones(length(tSvariants), 1));
            variantIndex = variantMatch_logic * transpose(1:length(tSvariants));
            variantIndex = variantIndex(isTwinning);
        end
        
        function [variantIndex] = variantDetermination2(tS, misOrientation, deviationInDegree)
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
        
        function [anglesMisoriFromTwin] = anglefromTwin(tS, misOrientation)
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
            anglesMisoriFromTwin = min(variantMatch_angle, [], 2);
        end

        function isC = isCompound(tS)
            % A twin is compound if the plane of shear S is a mirror plane.
            % The plane of shear S is the cross product of eta1 and k1.
            S = cross(tS.eta1, tS.k1);
            isC = S.isMirrorPlane;
        end

        function misori = parentTwinMisorientation(tS)
            % The misorientation depends on the twin type.
            misori = rotation.nan(size(tS));

            % Type I (1) and Compound (0)
            isType1OrCompound = tS.twinType == 1 | tS.twinType == 0;
            if any(isType1OrCompound)
                % Type I twin is defined by a reflection in the K1 plane.
                % For centrosymmetric crystals, this is equivalent to a 180-degree rotation about the normal to K1.
                misori(isType1OrCompound) = rotation.byAxisAngle(tS.k1(isType1OrCompound), pi);
            end
            % Type II (2)
            isType2 = tS.twinType == 2;
            if any(isType2)
                % Type II twin is defined by a 180-degree rotation about eta1.
                misori(isType2) = rotation.byAxisAngle(tS.eta1(isType2), pi);
            end
        end

        function shear = shearMagnitude(tS)
            shear = sqrt(4 * ((dot(tS.eta2.normalize, tS.k1.normalize).^-2) - 1));
        end
        
        function [shearvector, shearMagnitude] = shear(tS)
            shearvector = 2 * (tS.k1.normalize - dot(tS.eta2.normalize, tS.k1.normalize).^-1 * tS.eta2.normalize);
            shearvector.dispStyle = 'UVTW';
            shearMagnitude = norm(shearvector);
        end
        
        function defTensor = displacementGradient(tS)
            defTensor = (tS.shearMagnitude .* dyad(tS.eta1, tS.k1).normalize);
        end
        
        function correspondanceMatrix = correspondanceMatrix(tS)
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
                    dataCell = num2cell([d, reshape(tS.CRSS, [], 1)]);
                    fullRow = [dataCell, cellstr(typeNames)];
                    fullRow = fullRow.'; % Transpose for column-major printing
                    fprintf('eta_1 U   V   T   W | K_1 H   K   I   L   CRSS   Type\n');
                    fprintf('%7.0f %3.0f %3.0f %3.0f |  %3.0f %3.0f %3.0f %3.0f %6.0f   %s\n', fullRow{:});
                    % cprintf([d, reshape(tS.CRSS, [], 1), char({tS.twinType})], '-L', '  ', '-Lc', {'eta_1 U' 'V' 'T' 'W' '| K_1 H' 'K' 'I' 'L' 'CRSS' 'Type'});
                else
                    d = [tS.eta1.uvw tS.k1.hkl];
                    d(abs(d) < 1e-10) = 0; 
                    numericData = [d, reshape(tS.CRSS, [], 1)];
                    dataCell = [num2cell(numericData), cellstr(typeNames)];
                    cprintf(dataCell, '-L', '  ', '-Lc', {'eta_1 u' 'v' 'w' '| K_1 h' 'k' 'l' 'CRSS' 'Type'});
                    % cprintf([d, reshape(tS.CRSS, [], 1), char({tS.twinType})], '-L', '  ', '-Lc', {'eta_1 u' 'v' 'w' '| K_1 h' 'k' 'l' 'CRSS' 'Type'});
                end
            else
                d = round(100 * [tS.eta1.xyz tS.k1.xyz]) ./ 100;
                d(abs(d) < 1e-10) = 0;
                cprintf(d, '-L', '  ', '-Lc', {'x' 'y' 'z' ' |   x' 'y' 'z'});
            end
        end
        function nTS = mtimes(tS1, tS2)
            % Overloads the * operator to create a nestedTwinSystem
            % This allows for syntax like: doubleTwin = twin1 * twin2;
            nTS = nestedTwinSystem(tS1, tS2);
        end

        function n = length(tS)
            n = length(tS.k1);
        end

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