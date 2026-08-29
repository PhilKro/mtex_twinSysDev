classdef twinSystem
    properties
        CRSS     % Critical resolved shear stress
        eta1     % Shear direction
        eta2     % Conjugate shear direction
        k1       % Twin plane
        k2       % Conjugate twin plane
        rotAxis  % Zone axis
        twinType % Type of twin: 0:'Compound', 1:'Type I' or 2:'Type II'
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

                % Automatically determine twin type based on rational indices
                if nargin < 7 || isempty(type)
                    rK1 = twinSystem.isRational(tS.k1);
                    rK2 = twinSystem.isRational(tS.k2);

                    tS.twinType = ones(size(tS.k1)); % Default to Type I (1)
                    tS.twinType(rK2 & ~rK1) = 2;     % Type II
                    tS.twinType(rK1 & rK2) = 0;      % Compound
                    if any(~rK1 & ~rK2)
                        warning('For some twins both K1 and K2 plane were determined as irrational, thus no twin type could be determined. Defaulted to Type I.')
                    end
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

        function isC = isCompound(tS)
            % % A twin is compound if the plane of shear S is a mirror plane.
            % % The plane of shear S is the cross product of eta1 and k1.
            % S  = cross(tS.eta1, tS.k1);
            % isC = S.isMirrorPlane;
            isC = eq(tS.twinType, 0);
        end

        function misori = parentTwinMisorientation(tS)
            % The misorientation depends on the twin type.
            misori = orientation.nan(tS.CS, tS.CS, size(tS));

            % Type I (1) and Compound (0)
            isType1OrCompound = tS.twinType == 1 | tS.twinType == 0;
            if any(isType1OrCompound)
                % Type I twin is defined by a reflection in the K1 plane.
                % For centrosymmetric crystals, this is equivalent to a 180-degree rotation about the normal to K1.
                warning('twinSystem:parentTwinMisorientation:type1Assumption', ...
                    'The current implementation works with an improper rotation here. For orientation analysis this might not work (especially if you work with a centrosymmteric point group). Maybe the proper rotation should be implemented here aswell/instead?');
                misori(isType1OrCompound) = -orientation.byAxisAngle(tS.k1(isType1OrCompound), pi);
            end
            % Type II (2)
            isType2 = tS.twinType == 2;
            if any(isType2)
                % Type II twin is defined by a 180-degree rotation about eta1.
                misori(isType2) = orientation.byAxisAngle(tS.eta1(isType2), pi);
            end
        end

        function sM = shearMagnitude(tS)
            [~, sM] = tS.shear();
        end

        function [shearvector, shearMagnitude] = shear(tS)
            shearvector = 2 * (tS.k1.normalize - dot(tS.eta2.normalize, tS.k1.normalize).^-1 .* tS.eta2.normalize);
            shearvector.dispStyle = 'UVTW';
            shearMagnitude = norm(shearvector);
        end

        function defTensor = displacementGradient(tS, varargin)
            % DISPLACEMENTGRADIENT Calculates the displacement gradient (deformation tensor) of the twin
            % strictly in the local crystal reference frame.
            %
            % Syntax
            %   defTensor = displacementGradient(tS)
            [~, shearMagnitude] = tS.shear;
            defTensor = shearMagnitude .* dyad(tS.eta1, tS.k1).normalize;
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

                [~,shearMag] = shear(tS);

                if tS.eta1.lattice.isTriHex

                    c = tS.CS.cAxis;
                    delta_c = tS.displacementGradient('crystalCoords') * c;
                    val = dot(c, delta_c);


                    twinMode = strings(length(tS),1);
                    twinMode(val > 1e-4) = "extension";
                    twinMode(val < -1e-4) = "compression";
                    twinMode(val >= -1e-4 & val <= 1e-4) = "";

                    reta1 = round(tS.eta1);
                    rk1 = round(tS.k1);
                    d = [reta1.UVTW rk1.hkil];
                    d(abs(d) < 1e-10) = 0;

                    numericData = [d, reshape(tS.CRSS, [], 1), reshape(tS.variantId, [], 1), shearMag];
                    dataCell = [num2cell(numericData), cellstr(typeNames), cellstr(twinMode)];

                    cprintf(dataCell, '-L', '  ', '-Lc', {'eta_1 U' 'V' 'T' 'W' '| K_1 H' 'K' 'I' 'L' 'CRSS' 'Variant' 'Shear' 'Type' 'Mode'});

                else
                    reta1 = round(tS.eta1);
                    rk1 = round(tS.k1);
                    d = [reta1.uvw rk1.hkl];
                    d(abs(d) < 1e-10) = 0;
                    numericData = [d, reshape(tS.CRSS, [], 1), reshape(tS.variantId, [], 1), shearMag];
                    dataCell = [num2cell(numericData), cellstr(typeNames)];
                    cprintf(dataCell, '-L', '  ', '-Lc', {'eta_1 u' 'v' 'w' '| K_1 h' 'k' 'l' 'CRSS' 'Variant' 'Shear' 'Type'});
                end
            else
                d = round(100 * [tS.eta1.xyz tS.k1.xyz]) ./ 100;
                d(abs(d) < 1e-10) = 0;
                cprintf(d, '-L', '  ', '-Lc', {'x' 'y' 'z' ' |   x' 'y' 'z'});
            end
        end

        function out = mtimes(A, B)
            %MTIMES implement multiplication for twinSystem
            %
            % Syntax
            %   misori = tS1 * tS2
            %   ori = ori * tS
            %   ori = grains * tS

            if isa(A, 'twinSystem') && isa(B, 'twinSystem')
                % Returns the compounded geometric misorientation
                out = A.parentTwinMisorientation * B.parentTwinMisorientation;

            elseif (isa(A, 'orientation') || isa(A, 'grain2d')) && isa(B, 'twinSystem')
                if isa(A, 'grain2d')
                    oriA = A.meanOrientation;
                else
                    oriA = A;
                end

                % Reshape oriA to match B dimensions if needed, or MTEX broadcast handles it
                out = oriA .* B.parentTwinMisorientation;

            else
                out = builtin('mtimes', A, B);
            end
        end

        function F_total = calcDeformationSequence(tS_sym, variantPath, parentOri)
            % CALCDEFORMATIONSEQUENCE Computes the total deformation gradient
            % across a sequence of twin variants, optionally rotated into the specimen frame.
            %
            % Syntax
            %   F_total = tS_sym.calcDeformationSequence(variantPath)
            %   F_total = tS_sym.calcDeformationSequence(variantPath, parentOri)
            %
            % Note: This method expects tS_sym to be the symmetrized variant array.

            if isempty(variantPath)
                F_total = eye(3);
                return;
            end

            F_total = eye(3);

            % Compute crystal-frame sequential deformation
            for p = 1:length(variantPath)
                vIdx = variantPath(p);
                % Inside class methods, array indexing bypasses custom subsref.
                % We MUST use select() to safely index the vectorized object.
                vTS = tS_sym.select(vIdx);

                H_tensor = vTS.displacementGradient;
                if isa(H_tensor, 'tensor') || isa(H_tensor, 'tensor3')
                    H_mat = H_tensor.M;
                else
                    H_mat = double(H_tensor);
                end

                F_step = eye(3) + H_mat;
                F_total = F_step * F_total;
            end

            % Transform tensor from crystal to specimen coordinates using the parent orientation
            if nargin > 2 && ~isempty(parentOri)
                % Tensor rotation: R * F * R^T
                R = parentOri.matrix;
                F_total = R * F_total * R';
            end
        end


        function n = length(tS)
            if isempty(tS)
                n = 0;
            else
                n = length(tS.k1);
            end
        end

        function n = numel(tS,varargin)
            if isempty(tS)
                n = 0;
            else
                n = numel(tS.k1,varargin{:});
            end
        end

        function n = numArgumentsFromSubscript(tS, s, builtinContext)
            % Because tS is a vectorized object, we only ever want to return
            % a single object/array from any indexing operation (like tS.parent).
            % This prevents MATLAB from crashing with 'too few outputs' when numel > 1.
            n = 1;
        end

        function tS_out = select(tS, idx)
            % SELECT Helper method for safe array indexing inside class methods.
            % Bypasses native MATLAB bounds checking which treats tS as 1x1.
            s = struct('type', '()', 'subs', {{idx}});
            tS_out = subsref(tS, s);
        end

        function varargout = size(tS,varargin)
        % overloads size
            if builtin('numel', tS) == 0
                [varargout{1:nargout}] = builtin('size', tS, varargin{:});
            else
                [varargout{1:nargout}] = size(tS.eta1.x, varargin{:});
            end
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
            rotAxis = cross(k1, eta2);
            eta1 = cross(k1, rotAxis);
            k2 = cross(eta2, rotAxis);
            if nargin <3
                tS = twinSystem(rotAxis, k1, eta1, k2, eta2, 1);
            else
                tS = twinSystem(rotAxis, k1, eta1, k2, eta2, 1, type);
            end
            shearvec = 2 * (tS.k1.normalize - dot(tS.eta2.normalize, tS.k1.normalize).^-1 * tS.eta2.normalize);
            assert(dot(shearvec, tS.eta1, 'noSymmetry') > 0, "Something went wrong in the definition of the twin system: eta1 has the wrong sense. Sense of shear will not be correct")
        end

        function tS = byK2eta1(k2, eta1, type)
            if nargin < 3, type = []; end
            latticeIsHex = eq(eta1.CS.lattice, latticeType.hexagonal);

            rotAxis = cross(eta1, k2);
            if latticeIsHex
                rotAxis.dispStyle = 'UVTW';
            else
                rotAxis.dispStyle = 'uvw';
            end

            k1 = cross(rotAxis, eta1);
            if latticeIsHex
                k1.dispStyle = 'hkil';
            else
                k1.dispStyle = 'hkl';
            end

            eta2 = cross(rotAxis, k2);
            if latticeIsHex
                eta2.dispStyle = 'UVTW';
            else
                eta2.dispStyle = 'uvw';
            end

            tS = twinSystem(rotAxis, k1, eta1, k2, eta2, 1, type);

            shearvec = 2 * (tS.k1.normalize - dot(tS.eta2.normalize, tS.k1.normalize).^-1 * tS.eta2.normalize);
            assert(dot(shearvec, tS.eta1, 'noSymmetry') > 0, "Something went wrong in the definition of the twin system: eta1 has the wrong sense. Sense of shear will not be correct")
        end

        % function tS = byK1rotAxis(k1, rotAxis, type)
        %     error('twinSystem:byK1rotAxis:deprecated',"This method is not fully implemented yet. Use byK1eta2 or byK2eta1 instead for now.")
        %     if nargin < 3, type = 2; end
        %     if dot(k1, rotAxis, 'noSymmetry') < 1.00e-012
        %         eta1 = cross(k1, rotAxis);
        %         if eq(eta1.CS.lattice, latticeType.hexagonal)
        %             eta1.dispStyle = 'UVTW';
        %         else
        %             eta1.dispStyle = 'uvw';
        %         end
        %         k2 = eta1;
        %         if eq(eta1.CS.lattice, latticeType.hexagonal)
        %             k2.dispStyle = 'hkil';
        %         else
        %             k2.dispStyle = 'hkl';
        %         end
        %         eta2 = k1;
        %         eta2.dispStyle = 'UVTW';
        %         disp(append(newline, 'Chose a set of values for K2 and eta2!'))
        %         tS = twinSystem(rotAxis, k1, eta1, k2, eta2, 1, type);
        %     else
        %         error("k1 and rotAxis have to be orthogonal")
        %     end
        % end

        % function tS = hexagonal_110K(K, CS)
        %     if CS.lattice.isTriHex
        %         tS = arrayfun(@(x_int) twinSystem.byK1rotAxis(Miller({1, 0, -1, x_int}, CS), Miller({1, -2, 1, 0}, CS, 'UVTW'), 2), K);
        %     else
        %         disp("Provide a hexagonal CS for this method")
        %     end
        % end
        %
        % function tS = hexagonal_211K(K, CS)
        %     if CS.lattice.isTriHex
        %         tS = arrayfun(@(x_int) twinSystem.byK1rotAxis(Miller({1, 1, -2, x_int}, CS), Miller({1, -1, 0, 0}, CS, 'UVTW'), 2), K);
        %     else
        %         disp("Provide a hexagonal CS for this method")
        %     end
        % end

        function tS = hexagonal_1012(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.calculateTheoreticalTwins(CS,Miller(1,0,-1,2,CS),2);
                if ~isempty(tS), idx = find(tS.twinType == 0, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_2241(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.calculateTheoreticalTwins(CS,Miller(2,2,-4,1,CS),2);
                if ~isempty(tS), idx = find(tS.twinType == 0, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1011(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.calculateTheoreticalTwins(CS,Miller(1,0,-1,1,CS),4);
                if ~isempty(tS), idx = find(tS.twinType == 0, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1011_TypeI(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.calculateTheoreticalTwins(CS,Miller(1,0,-1,1,CS),2);
                if ~isempty(tS), idx = find(tS.twinType == 1, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_2021(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.calculateTheoreticalTwins(CS,Miller(2,0,-2,1,CS),2);
                if ~isempty(tS), idx = find(tS.twinType == 0, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1121(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.calculateTheoreticalTwins(CS,Miller(1,1,-2,1,CS),1);
                if ~isempty(tS), idx = find(tS.twinType == 0, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1013(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.calculateTheoreticalTwins(CS,Miller(1,0,-1,3,CS),4);
                if ~isempty(tS), idx = find(tS.twinType == 0, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1013_TypeI(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.calculateTheoreticalTwins(CS,Miller(1,0,-1,3,CS),2);
                if ~isempty(tS), idx = find(tS.twinType == 1, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1122(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.calculateTheoreticalTwins(CS,Miller(1,1,-2,2,CS),3);
                if ~isempty(tS), idx = find(tS.twinType == 0, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1124(CS)
            if CS.lattice.isTriHex
                tS = twinSystem.calculateTheoreticalTwins(CS,Miller(1,1,-2,4,CS),3);
                if ~isempty(tS), idx = find(tS.twinType == 0, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = bct_101(CS)
            % classic fcc {111}<112> twin system
            if eq(CS.lattice, latticeType.tetragonal)
                tS = twinSystem.calculateTheoreticalTwins(CS, Miller(1,0,1,CS,'hkl'), 4, 'BCC');
                if ~isempty(tS), idx = find(tS.twinType == 0, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a tetragonal CS for this method")
            end
        end

        function tS = fcc_111(CS)
            % classic fcc {111}<112> twin system
            if eq(CS.lattice, latticeType.cubic)
                tS = twinSystem.calculateTheoreticalTwins(CS, Miller(1,1,1,CS,'hkl'), 2, 'FCC');
                if ~isempty(tS), idx = find(tS.twinType == 0, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a cubic CS for this method")
            end
        end

        function tS = bcc_112(CS)
            % classic bcc {112}<111> twin system
            if eq(CS.lattice, latticeType.cubic)
                tS = twinSystem.calculateTheoreticalTwins(CS, Miller(1,1,2,CS,'hkl'), 2, 'BCC');
                if ~isempty(tS), idx = find(tS.twinType == 0, 1); if ~isempty(idx), tS = tS.select(idx); else, tS = twinSystem.empty; end; end
            else
                disp("Provide a cubic CS for this method")
            end
        end
    function printTheoreticalTwins(twins)
            % PRINTTHEORETICALTWINS Prints an array of twin structures as a formatted table.
            %
            % theotw = output struct from calculateTheoreticalTwins
            % Syntax
            %   twinSystem.printTheoreticalTwins(theotw)

            if isempty(twins)
                disp('No theoretical twins found.');
                return;
            end

            % Define column widths
            wNo = 5; wK1 = 15; wK2 = 15; wEta1 = 15; wEta2 = 15; wQ = 4; wShear = 10;

            % Print Header
            fprintf('\n%-*s | %-*s | %-*s | %-*s | %-*s | %-*s | %-*s\n', ...
                wNo, 'No.', wK1, 'K1 (Plane)', wK2, 'K2 (Plane)', ...
                wEta1, 'eta1 (Dir)', wEta2, 'eta2 (Dir)', wQ, 'q', wShear, 'Shear');

            % Print Separator
            totalWidth = wNo + wK1 + wK2 + wEta1 + wEta2 + wQ + wShear + 18;
            fprintf('%s\n', repmat('-', 1, totalWidth));

            % Print Rows
            for i = 1:length(twins)
                % Extract and format Miller indices
                % char() converts the MTEX object to a string, strtrim removes trailing spaces
                dispK1 = twins(i).K1; dispK1.dispStyle = 'hkl';
                dispK2 = twins(i).K2; dispK2.dispStyle = 'hkl';
                dispEta1 = twins(i).eta1; dispEta1.dispStyle = 'uvw';
                dispEta2 = twins(i).eta2; dispEta2.dispStyle = 'uvw';

                k1_str   = strtrim(char(round(dispK1,'maxHKL',100),'hkl'));
                k2_str   = strtrim(char(round(dispK2,'maxHKL',100),'hkl'));
                eta1_str = strtrim(char(round(dispEta1,'maxHKL',100),'uvw'));
                eta2_str = strtrim(char(round(dispEta2,'maxHKL',100),'uvw'));

                % Optionally wrap in standard crystallographic brackets if MTEX char() doesn't
                if ~startsWith(k1_str, '(') && ~startsWith(k1_str, '{')
                    k1_str = sprintf('(%s)', k1_str);
                    k2_str = sprintf('(%s)', k2_str);
                    eta1_str = sprintf('[%s]', eta1_str);
                    eta2_str = sprintf('[%s]', eta2_str);
                end

                % Print formatted row
                fprintf('%-*d | %-*s | %-*s | %-*s | %-*s | %-*d | %-*.4f\n', ...
                    wNo, i, wK1, k1_str, wK2, k2_str, ...
                    wEta1, eta1_str, wEta2, eta2_str, wQ, twins(i).q, wShear, twins(i).shear);
            end
            fprintf('\n');
        end

        function [twin_objects, twin_structs] = calculateTheoreticalTwins(cs, inputObjOrMaxIndex, qRange, structName, varargin)
            %CALCULATETHEORETICALTWINS Generates possible twin systems
            %(normal mode) for a given crystalSymmetry, plane or max hkl
            %and twin index (q). 'FCC' or 'BCC' can be handed over as char.
            %
            % Acknowledgement:
            %   The logic and algorithms in this method are translated and adapted
            %   from the ARPGE, GenOVa, and Crystals programs developed by Cyril Cayron.
            %
            % Inputs
            %   cs                  - crystalSymmetry object defining the lattice symmetry.
            %   inputObjOrMaxIndex  - Miller object representing the candidate twin plane (K1), direction (eta1), or maximum Miller index.
            %   qRange              - Double array defining the range of twin indices to evaluate.
            %   structName          - String defining the lattice type (e.g., 'FCC', 'BCC', 'none') for centering rules.
            %   maxShear            - (Optional) Maximum allowed shear magnitude (default: 1.0).
            %   twinTypes           - (Optional) Char/Cell array specifying what to calculate: 'type1', 'type2', or both.
            %
            % Outputs
            %   twin_objects        - Array of initialized twinSystem objects for valid modes.
            %   twin_structs        - Structured array containing the detailed crystallographic
            %                           elements (K1, K2, eta1, eta2, shear, q, F).

            if nargin < 4 || isempty(structName), structName = 'none'; end % Default structure
            if nargin < 3 || isempty(qRange), qRange = 1:4; end % Default qRange

            % Parse optional varargin for maxShear and twinTypes
            maxShear = 1.0;
            twinTypes = {};
            if nargin >= 5
                if isnumeric(varargin{1})
                    maxShear = varargin{1};
                    if length(varargin) > 1
                        twinTypes = varargin{2};
                    end
                else
                    twinTypes = varargin{1};
                    if length(varargin) > 1 && isnumeric(varargin{2})
                        maxShear = varargin{2};
                    end
                end
            end

            % If twinTypes is provided as char/string, convert to cell
            if ischar(twinTypes) || isstring(twinTypes)
                twinTypes = cellstr(twinTypes);
            end

            twin_objects = twinSystem.empty;
            twin_structs = struct.empty;

            if isa(inputObjOrMaxIndex, 'Miller')
                inputObj = inputObjOrMaxIndex;
                if inputObj.CS ~= cs
                    error('Crystal symmetry of the input object must match the provided crystal symmetry.');
                end

                isReciprocal = strcmpi(inputObj.dispStyle, 'hkl') || strcmpi(inputObj.dispStyle, 'hkil');
                isDirect = strcmpi(inputObj.dispStyle, 'uvw') || strcmpi(inputObj.dispStyle, 'UVTW');

                % If twinTypes is specified, check against what is physically possible
                if ~isempty(twinTypes)
                    forceType1 = any(strcmpi('type1', twinTypes) | strcmpi('Type1', twinTypes));
                    forceType2 = any(strcmpi('type2', twinTypes) | strcmpi('Type2', twinTypes));

                    if isReciprocal && forceType2
                        warning('Cannot calculate Type 2 twins for an input plane. Only Type 1 will be calculated.');
                    elseif isDirect && forceType1
                        warning('Cannot calculate Type 1 twins for an input direction. Only Type 2 will be calculated.');
                    end
                end

                if isReciprocal
                    [twin_objects, twin_structs] = twinSystem.calculateForPlane(inputObj, qRange, structName, maxShear);
                elseif isDirect
                    [twin_objects, twin_structs] = twinSystem.calculateForDirection(inputObj, qRange, structName, maxShear);
                else
                    error('Unrecognized Miller object type.');
                end

            elseif isnumeric(inputObjOrMaxIndex)
                maxMillerIndex = inputObjOrMaxIndex;

                % Default to both types if a max index is passed and twinTypes isn't specified
                if isempty(twinTypes)
                    twinTypes = {'type1', 'type2'};
                end

                calcType1 = any(strcmpi('type1', twinTypes) | strcmpi('Type1', twinTypes));
                calcType2 = any(strcmpi('type2', twinTypes) | strcmpi('Type2', twinTypes));

                % Generate grid of indices
                [H, K, L] = ndgrid(-maxMillerIndex:maxMillerIndex, -maxMillerIndex:maxMillerIndex, -maxMillerIndex:maxMillerIndex);
                indices = [H(:), K(:), L(:)];
                indices(all(indices == 0, 2), :) = [];

                for i = 1:size(indices, 1)
                    idx = indices(i, :);

                    if calcType1
                        plane = Miller(idx(1), idx(2), idx(3), cs, 'hkl');
                        [plane_twin_objects, plane_twin_structs] = twinSystem.calculateForPlane(plane, qRange, structName, maxShear);
                        if ~isempty(plane_twin_objects)
                            twin_objects = [twin_objects; plane_twin_objects]; %#ok<AGROW>
                            twin_structs = [twin_structs; plane_twin_structs]; %#ok<AGROW>
                        end
                    end

                    if calcType2
                        dir = Miller(idx(1), idx(2), idx(3), cs, 'uvw');
                        [dir_twin_objects, dir_twin_structs] = twinSystem.calculateForDirection(dir, qRange, structName, maxShear);
                        if ~isempty(dir_twin_objects)
                            twin_objects = [twin_objects; dir_twin_objects]; %#ok<AGROW>
                            twin_structs = [twin_structs; dir_twin_structs]; %#ok<AGROW>
                        end
                    end
                end
            else
                error('Second argument must be a Miller object (for a single plane/direction) or a numeric value (for max Miller index).');
            end
        end

       function plotShearPlane(plane, qRange, structName, defaultQ, plotTwinned)
            % PLOTSHEARPLANE Visualizes the 2D plane of shear for a given twin system.
            %
            % This method generates a 2D projection of the crystal lattice strictly
            % along the plane of shear defined by a candidate twin plane (K1) and its
            % calculated shear direction (eta1). By optionally plotting both the
            % parent and twinned lattices, it can generate a 2D dichromatic pattern.
            %
            % The resulting figure includes an interactive dropdown menu allowing you
            % to switch the projection between different theoretical twin modes found.

            % Ensure default variables are set
            if nargin < 5, plotTwinned = false; end
            if nargin < 4, defaultQ = 1; end
            if nargin < 3, structName = 'none'; end
            if nargin < 2, qRange = 1:4; end

            cs = plane.CS;
            isHCP = strcmpi(structName, 'HCP') && cs.lattice.isTriHex;

            [~, twin_structs] = twinSystem.calculateForPlane(plane, qRange, structName);
            maxQ = max(qRange);

            if isempty(twin_structs)
                warning('No theoretical twins found to plot.');
                return;
            end

            h = round(plane.h); k = round(plane.k); l = round(plane.l);
            [~, G_rec] = twinSystem.getMetricTensors(cs);

            [OZ1, U, V] = twinSystem.Bezout3D(plane);
            [~, U, V] = twinSystem.findUV(0.5*U, 0.5*V, structName);

            ap = [h; k; l];
            apD = G_rec * ap;
            dhkl2 = 1 / (ap' * apD);
            OH1_vec = dhkl2 * apD;

            [G_dir, ~] = twinSystem.getMetricTensors(cs);

            % Generate the 3D supercell once
            grid_range = 60;
            [U_grid, V_grid, M_grid] = ndgrid(-grid_range:grid_range, -grid_range:grid_range, 0:maxQ+2);
            atoms3D = U.uvw' * U_grid(:)' + V.uvw' * V_grid(:)' + OZ1.uvw' * M_grid(:)';

            if isHCP
                hcp_shift = [1/3; 2/3; 1/2];
                hcp_atoms3D = atoms3D + hcp_shift;
            end

            % Prepare UI labels for the modes
            mode_labels = cell(1, length(twin_structs));
            for i = 1:length(twin_structs)
                k2_str = char(round(twin_structs(i).K2, 'maxHKL', 100));
                eta2_str = char(round(twin_structs(i).eta2, 'maxHKL', 100));
                % Strip any existing brackets from MTEX char output
                k2_str = strtrim(strrep(strrep(strrep(strrep(k2_str, '(', ''), ')', ''), '{', ''), '}', ''));
                eta2_str = strtrim(strrep(strrep(strrep(strrep(eta2_str, '[', ''), ']', ''), '<', ''), '>', ''));

                mode_labels{i} = sprintf('q=%d, S=%.3f, K2=(%s), eta2=[%s]', twin_structs(i).q, twin_structs(i).shear, k2_str, eta2_str);
            end

            % Find initial index
            idx = find([twin_structs.q] == defaultQ, 1);
            if isempty(idx)
                warning('Default q=%d not found in evaluated modes. Using q=%d as default.', defaultQ, twin_structs(1).q);
                idx = 1;
            end

            % Setup Figure and UI
            fig = figure('Name', 'Interactive Plane of Shear', 'Units', 'normalized', 'Position', [0.2 0.2 0.6 0.6]);

            % Layout: Top panel for UI, bottom axes
            top_panel = uipanel(fig, 'Position', [0 0.9 1 0.1], 'BorderType', 'none');
            ax = axes(fig, 'Position', [0.05 0.05 0.9 0.8]);

            uicontrol(top_panel, 'Style', 'text', 'String', 'Select Plane of Shear Mode:', ...
                'Units', 'normalized', 'Position', [0.02 0.3 0.25 0.4], 'HorizontalAlignment', 'right', 'FontSize', 11);

            dropdown = uicontrol(top_panel, 'Style', 'popupmenu', 'String', mode_labels, ...
                'Units', 'normalized', 'Position', [0.28 0.3 0.5 0.4], 'Value', idx, 'FontSize', 11, ...
                'Callback', @(~,~) drawPlot());

            % Draw initial plot
            drawPlot();

            function drawPlot()
                target_idx = dropdown.Value;
                target_twin = twin_structs(target_idx);

                eta1_vec = target_twin.eta1.uvw';
                norm_eta1 = sqrt(eta1_vec' * G_dir * eta1_vec);
                projX = @(vec) (vec' * G_dir * eta1_vec) / norm_eta1;

                norm_OH1 = sqrt(OH1_vec' * G_dir * OH1_vec);
                projY = @(vec) (vec' * G_dir * OH1_vec) / norm_OH1;

                plane_of_shear_normal = cross(eta1_vec, OH1_vec);
                pos_normal_norm = norm(plane_of_shear_normal);

                % Filter atoms exactly on the plane of shear
                dists = abs(plane_of_shear_normal' * atoms3D) / pos_normal_norm;
                valid_idx = dists < 1e-4;
                valid_atoms = atoms3D(:, valid_idx);

                X_all = (eta1_vec' * G_dir * valid_atoms) / norm_eta1;
                Y_all = (OH1_vec' * G_dir * valid_atoms) / norm_OH1;

                if isHCP
                    dists_hcp = abs(plane_of_shear_normal' * hcp_atoms3D) / pos_normal_norm;
                    valid_hcp_idx = dists_hcp < 1e-4;
                    valid_hcp_atoms = hcp_atoms3D(:, valid_hcp_idx);

                    X_hcp_all = (eta1_vec' * G_dir * valid_hcp_atoms) / norm_eta1;
                    Y_hcp_all = (OH1_vec' * G_dir * valid_hcp_atoms) / norm_OH1;
                end

                if plotTwinned
                    twinned_atoms3D = target_twin.F * valid_atoms;
                    X_tw_all = (eta1_vec' * G_dir * twinned_atoms3D) / norm_eta1;
                    Y_tw_all = (OH1_vec' * G_dir * twinned_atoms3D) / norm_OH1;

                    if isHCP
                        twinned_hcp_atoms3D = target_twin.F * valid_hcp_atoms;
                        X_hcp_tw_all = (eta1_vec' * G_dir * twinned_hcp_atoms3D) / norm_eta1;
                        Y_hcp_tw_all = (OH1_vec' * G_dir * twinned_hcp_atoms3D) / norm_OH1;
                    end
                end

                total_y_height = projY(maxQ * OH1_vec);

                all_x_ends = zeros(1, length(twin_structs));
                for i = 1:length(twin_structs)
                    all_x_ends(i) = projX(twin_structs(i).eta2.uvw');
                end

                X_min_base = min([0, all_x_ends]);
                X_max_base = max([0, all_x_ends]);

                padding_x = max(1.5 * total_y_height, 2 * norm_eta1);
                X_min = X_min_base - padding_x;
                X_max = X_max_base + padding_x;

                plot_idx = (X_all >= X_min) & (X_all <= X_max);
                X_atoms = X_all(plot_idx);
                Y_atoms = Y_all(plot_idx);

                if isHCP
                    plot_hcp_idx = (X_hcp_all >= X_min) & (X_hcp_all <= X_max);
                    X_hcp = X_hcp_all(plot_hcp_idx);
                    Y_hcp = Y_hcp_all(plot_hcp_idx);
                end

                if plotTwinned
                    plot_tw_idx = (X_tw_all >= X_min) & (X_tw_all <= X_max);
                    X_tw_atoms = X_tw_all(plot_tw_idx);
                    Y_tw_atoms = Y_tw_all(plot_tw_idx);

                    if isHCP
                        plot_hcp_tw_idx = (X_hcp_tw_all >= X_min) & (X_hcp_tw_all <= X_max);
                        X_hcp_tw = X_hcp_tw_all(plot_hcp_tw_idx);
                        Y_hcp_tw = Y_hcp_tw_all(plot_hcp_tw_idx);
                    end
                end

                cla(ax);
                hold(ax, 'on');

                scatter(ax, X_atoms, Y_atoms, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none');
                if isHCP
                    scatter(ax, X_hcp, Y_hcp, '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none');
                end

                if plotTwinned
                    scatter(ax, X_tw_atoms, Y_tw_atoms, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
                    if isHCP
                        scatter(ax, X_hcp_tw, Y_hcp_tw, '^', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
                    end
                end

                plot(ax, [X_min, X_max], [0, 0], 'r-', 'LineWidth', 2);

                for q = 1:maxQ
                    y_height = projY(q * OH1_vec);
                    plot(ax, [X_min, X_max], [y_height, y_height], 'r--');
                end

                max_y_height = projY(maxQ * OH1_vec);
                plot(ax, [0, 0], [0, max_y_height], 'k-');

                axis(ax, 'equal');
                title(ax, sprintf('Plane of Shear Projection: K1 = %s', char(plane)));
                xlabel(ax, sprintf('Shear Direction (\\eta_1) aligned to %s mode', mode_labels{target_idx}));
                ylabel(ax, 'Plane Normal (OH_1)');

                for i = 1:length(twin_structs)
                    OA = twin_structs(i).eta2.uvw';
                    x_end = projX(OA);
                    y_end = projY(OA);

                    current_eta1 = twin_structs(i).eta1.uvw';
                    current_pos_normal = cross(current_eta1, OH1_vec);

                    n1 = plane_of_shear_normal / pos_normal_norm;
                    n2 = current_pos_normal / norm(current_pos_normal);

                    try
                        k2_tex = twin_structs(i).K2.latexChar('round');
                        eta2_tex = twin_structs(i).eta2.latexChar('round');
                    catch
                        k2_tex = ['(' char(round(twin_structs(i).K2, 'maxHKL', 100)) ')'];
                        eta2_tex = ['[' char(round(twin_structs(i).eta2, 'maxHKL', 100)) ']'];
                    end

                    if abs(abs(dot(n1, n2)) - 1) < 1e-4
                        col = 'b';
                        ls = '-';
                        lbl_tex = sprintf(' $q=%d$, $S: %.3f$, $K_2: %s$, $\\eta_2: %s$', twin_structs(i).q, twin_structs(i).shear, k2_tex, eta2_tex);
                        valign = 'bottom';
                    else
                        col = 'r';
                        ls = '--';
                        lbl_tex = sprintf(' $q=%d$ (Diff PoS)', twin_structs(i).q);
                        valign = 'top';
                    end

                    plot(ax, [0, x_end], [0, y_end], 'Color', col, 'LineStyle', ls, 'LineWidth', 1.5);

                    vec_len = sqrt(x_end^2 + y_end^2);
                    if vec_len > 1e-4
                        head_L = min((X_max - X_min) * 0.02, vec_len * 0.3);
                        head_W = head_L * 0.75;

                        theta_arrow = atan2(y_end, x_end);
                        P_base = [0, 0; -head_L, head_W/2; -head_L, -head_W/2];
                        R = [cos(theta_arrow), -sin(theta_arrow); sin(theta_arrow), cos(theta_arrow)];
                        P_rot = (R * P_base')';
                        P_rot(:,1) = P_rot(:,1) + x_end;
                        P_rot(:,2) = P_rot(:,2) + y_end;

                        fill(ax, P_rot(:,1), P_rot(:,2), col, 'EdgeColor', 'none');
                    end

                    text(ax, x_end, y_end, lbl_tex, 'VerticalAlignment', valign, 'Color', col, 'Interpreter', 'latex', 'FontSize', 11);
                end

                xlim(ax, [X_min, X_max]);
                ylim(ax, [-0.5 * norm_OH1, max_y_height + norm_OH1]);

                hold(ax, 'off');
            end
        end
    end

    methods (Static, Access = private)

        %% --- UTILITY MATH FUNCTIONS ---

        function res = is_int(x)
            % Check if number(s) is an integer within floating point tolerance
            res = all(abs(x - round(x)) < 1e-6);
        end

        function res = is_sameodd(p)
            % Check if the coordinates are all odd, or all even
            d = 0.5;
            ceven = all(twinSystem.is_int(p * d));
            codd = all(twinSystem.is_int((p - 1) * d));
            res = ceven || codd;
        end

        %% --- SELECTION RULES ---

        function res = is_FCC_Dir(dir)
            uvw = dir.uvw;
            res = twinSystem.is_int(sum(uvw));
        end

        function res = is_BCC_Dir(dir)
            uvw = dir.uvw;
            res = twinSystem.is_sameodd(2 * uvw);
        end

        function res = is_FCC_Rec(plane)
            hkl = plane.hkl;
            res = twinSystem.is_sameodd(hkl);
        end

        function res = is_BCC_Rec(plane)
            hkl = plane.hkl;
            res = twinSystem.is_int(sum(hkl) * 0.5);
        end

        %% --- DIOPHANTINE SOLVER ---

        function [OZ1, U, V] = Bezout3D(plane)
            % BEZOUT3D Solves the 3D linear Diophantine equation for a crystallographic plane.
            %
            % Syntax
            %   [OZ1, U, V] = twinSystem.Bezout3D(plane)
            %
            % Description
            %   Applies the Extended Euclidean Algorithm to a plane's Miller indices (hkl)
            %   to calculate three physical, integer-based direct lattice vectors.
            %
            % Outputs
            %   OZ1 - Direct lattice vector connecting the origin to an atom on the
            %         first adjacent (hkl) plane.
            %   U   - First primitive lattice vector lying exactly within the (hkl) plane.
            %   V   - Second primitive lattice vector lying exactly within the (hkl) plane,
            %         linearly independent from U.
            hkl = round(plane.hkl);
            h = hkl(1); k = hkl(2); l = hkl(3);
            CS = plane.CS;

            % Extended Euclidean Algorithm
            [g, a, b] = gcd(h, k);
            [G, c, d] = gcd(g, l);

            % OZ1 is a solution such that h*u + k*v + l*w = G (ideally 1 for primitive planes)
            % the euclidian extension simplifies to (c*a)*h + (c*b)*k + d*l = G
            % these values are calced by gcd and assigned to uvw
            u = c * a;
            v = c * b;
            w = d;

            % OZ1 is a lattice vector pointing to an atom on an adjacent
            % hkl plane
            OZ1 = Miller(u, v, w, CS, 'uvw');

            % U and V form a basis for the plane (U.hkl = 0, V.hkl = 0)
            if g ~= 0
                U = Miller(k/g, -h/g, 0, CS, 'uvw');
                V = Miller(a*l/G, b*l/G, -g/G, CS, 'uvw');
            else
                % If h=0 and k=0, then the plane is (0,0,l)
                U = Miller(1, 0, 0, CS, 'uvw');
                V = Miller(0, 1, 0, CS, 'uvw');
            end
        end

        function [OZ1, H, K] = Bezout3DRec(dir)
            % BEZOUT3DREC Solves the 3D linear Diophantine equation for a crystallographic direction.
            %
            % Syntax
            %   [OZ1, H, K] = twinSystem.Bezout3DRec(dir)
            %
            % Description
            %   Applies the Extended Euclidean Algorithm to a direction's Miller indices (uvw)
            %   to calculate three physical, integer-based reciprocal lattice vectors.
            %
            % Outputs
            %   OZ1 - Reciprocal lattice vector connecting the origin to an atom on the
            %         first adjacent layer perpendicular to the direction.
            %   H   - First primitive reciprocal lattice vector lying exactly perpendicular to dir.
            %   K   - Second primitive reciprocal lattice vector lying exactly perpendicular to dir,
            %         linearly independent from H.
            uvw = round(dir.uvw);
            u = uvw(1); v = uvw(2); w = uvw(3);
            CS = dir.CS;

            % Extended Euclidean Algorithm
            [g, a, b] = gcd(u, v);
            [G, c, d] = gcd(g, w);

            h = c * a;
            k = c * b;
            l = d;

            OZ1 = Miller(h, k, l, CS, 'hkl');

            if g ~= 0
                H = Miller(v/g, -u/g, 0, CS, 'hkl');
                K = Miller(a*w/G, b*w/G, -g/G, CS, 'hkl');
            else
                % If u=0 and v=0, then the direction is (0,0,w)
                H = Miller(1, 0, 0, CS, 'hkl');
                K = Miller(0, 1, 0, CS, 'hkl');
            end
        end

        %% --- BASIS VECTOR MINIMIZATION ---

        function [oz2, Ures, Vres] = findUV(U, V, structName)
            % Finds the smallest integer/half-integer directions belonging to the U,V plane
            CS = U.CS;
            oz2 = Miller(0, 0, 0, CS, 'uvw');

            if strcmpi(structName, 'FCC')
                u_fcc = twinSystem.is_FCC_Dir(U);
                v_fcc = twinSystem.is_FCC_Dir(V);
                if u_fcc && v_fcc
                    Ures = U; Vres = V;
                elseif u_fcc && ~v_fcc
                    Ures = U; Vres = 2 * V; oz2 = V;
                elseif v_fcc && ~u_fcc
                    Ures = 2 * U; Vres = V; oz2 = U;
                else
                    Ures = U - V; Vres = U + V; oz2 = U;
                end
            elseif strcmpi(structName, 'BCC')
                u_bcc = twinSystem.is_BCC_Dir(U);
                v_bcc = twinSystem.is_BCC_Dir(V);
                if u_bcc && v_bcc
                    Ures = U; Vres = V;
                elseif u_bcc && ~v_bcc
                    Ures = U; Vres = 2 * V; oz2 = V;
                elseif v_bcc && ~u_bcc
                    Ures = 2 * U; Vres = V; oz2 = U;
                else
                    Ures = U - V; Vres = U + V; oz2 = U;
                end
            else
                Ures = U; Vres = V;
            end
        end

        function [oz2, Hres, Kres] = findHK(H, K, structName)
            % Finds the smallest integer/half-integer planes belonging to the H,K zone
            CS = H.CS;
            oz2 = Miller(0, 0, 0, CS, 'hkl');

            if strcmpi(structName, 'FCC')
                h_fcc = twinSystem.is_FCC_Rec(H);
                k_fcc = twinSystem.is_FCC_Rec(K);
                if h_fcc && k_fcc
                    Hres = H; Kres = K;
                elseif h_fcc && ~k_fcc
                    Hres = H; Kres = 2 * K; oz2 = K;
                elseif k_fcc && ~h_fcc
                    Hres = 2 * H; Kres = K; oz2 = H;
                else
                    Hres = H - K; Kres = H + K; oz2 = H;
                end
            elseif strcmpi(structName, 'BCC')
                h_bcc = twinSystem.is_BCC_Rec(H);
                k_bcc = twinSystem.is_BCC_Rec(K);
                if h_bcc && k_bcc
                    Hres = H; Kres = K;
                elseif h_bcc && ~k_bcc
                    Hres = H; Kres = 2 * K; oz2 = K;
                elseif k_bcc && ~h_bcc
                    Hres = 2 * H; Kres = K; oz2 = H;
                else
                    Hres = H - K; Kres = H + K; oz2 = H;
                end
            else
                Hres = H; Kres = K;
            end
        end

        %% --- SHEAR HELPER METRICS ---

        function [G_dir, G_rec] = getMetricTensors(CS)
            % Builds the direct and reciprocal metric tensors dynamically from MTEX axes
            a = CS.aAxis.xyz';
            b = CS.bAxis.xyz';
            c = CS.cAxis.xyz';

            % 3x3 matrix where columns are orthonormal spatial basis vectors
            M = [a, b, c];
            G_dir = M' * M;
            G_rec = inv(G_dir);
        end

        function bool = isRational(m, maxRationalIndex)
            if nargin < 2, maxRationalIndex = 5; end
            [~, err] = round(m, 'maxHKL', 2 * maxRationalIndex, 'hkl');
            bool = err < 1e-4;
        end

        function res = ShearBevisCrocker(C, CS)
            % C = correspondence matrix
            % [G_dir, G_rec] = twinSystem.getMetricTensors(CS);
            [G_dir, G_rec] = twinSystem.getMetricTensors(CS);
            res2 = trace(G_rec * C' * G_dir * C) - 3;
            if res2 >= 0
                res = sqrt(res2);
            else
                res = 0;
            end
        end

        function res = ShearMyFormula(F, CS)
            % F = distortion matrix
            % [G_dir, G_rec] = twinSystem.getMetricTensors(CS);
            [G_dir, G_rec] = twinSystem.getMetricTensors(CS);
            I = eye(3);
            res = sqrt(trace(G_dir * (F - I) * G_rec * (F - I)'));
        end

        function res = ShearMyFormulaB(F, CS)
            % [G_dir, ~] = twinSystem.getMetricTensors(CS);
            [G_dir, ~] = twinSystem.getMetricTensors(CS);
            I = eye(3);
            res = sqrt(trace((F - I)' * G_dir * (F - I)));
        end

        function res = ShearMyFormulaC(F, CS)
            % [G_dir, G_rec] = twinSystem.getMetricTensors(CS);
            [G_dir, G_rec] = twinSystem.getMetricTensors(CS);
            res = sqrt(trace(G_dir * F * G_rec * F') - 3);
        end

        function [plane_twin_objects, plane_twin_structs] = calculateForPlane(plane, qRange, structName, maxShear)
            % CALCULATEFORPLANE Determines the lowest shear normal Type I (or compound) twinning modes for a given plane and range of twin indeces (qRange).
            %
            % Acknowledgement:
            %   The logic and algorithms in this method are translated and adapted
            %   from the ARPGE, GenOVa, and Crystals programs developed by Cyril Cayron.
            % Syntax
            %   [plane_twin_objects, plane_twin_structs] = twinSystem.calculateForPlane(plane, qRange, structName)
            %
            % Description
            %   Evaluates a candidate twin plane (K1) to find viable conjugate shear directions
            %   (eta2) up to a maximum twin index (q). By using the physical
            %   lattice node closest to the geometric normal on the q-th layer, this function
            %   isolates "normal modes" as defined by C. Cayron. https://www.mdpi.com/2075-4701/10/2/231
            %
            % Inputs
            %   plane      - Miller object representing the candidate twin plane (K1).
            %   qRange     - Double array defining the range of twin indices to evaluate.
            %   structName - String defining the lattice type (e.g., 'FCC', 'BCC', 'none') for centering rules.
            %
            % Outputs
            %   plane_twin_objects - Array of initialized twinSystem objects for valid modes.
            %   plane_twin_structs - Structured array containing the detailed crystallographic
            %                        elements (K1, K2, eta1, eta2, shear, q, F).

            if nargin < 4 || isempty(maxShear), maxShear = 1.0; end
            [plane_twin_objects, plane_twin_structs] = twinSystem.calculateTheoreticalMode(plane, qRange, structName, 'Type I', maxShear);
        end

        function [dir_twin_objects, dir_twin_structs] = calculateForDirection(dir, qRange, structName, maxShear)
            % CALCULATEFORDIRECTION Determines the lowest shear normal Type II (or compound) twinning modes for a given direction and range of twin indeces (qRange).
            %
            % Acknowledgement:
            %   The logic and algorithms in this method are translated and adapted
            %   from the ARPGE, GenOVa, and Crystals programs developed by Cyril Cayron.
            % Syntax
            %   [dir_twin_objects, dir_twin_structs] = twinSystem.calculateForDirection(dir, qRange, structName)
            %
            % Description
            %   Evaluates a candidate twin direction (eta1) to find viable conjugate twin planes
            %   (K2) up to a maximum twin index (q). By using the physical
            %   lattice node closest to the geometric normal on the q-th layer in reciprocal space,
            %   this function isolates "normal modes" as defined by C. Cayron.
            %   https://www.mdpi.com/2075-4701/10/2/231
            %
            % Inputs
            %   dir        - Miller object representing the candidate twin direction (eta1).
            %   qRange     - Double array defining the range of twin indices to evaluate.
            %   structName - String defining the lattice type (e.g., 'FCC', 'BCC', 'none') for centering rules.
            %
            % Outputs
            %   dir_twin_objects - Array of initialized twinSystem objects for valid modes.
            %   dir_twin_structs - Structured array containing the detailed crystallographic
            %                      elements (K1, K2, eta1, eta2, shear, q, F).

            if nargin < 4 || isempty(maxShear), maxShear = 1.0; end
            [dir_twin_objects, dir_twin_structs] = twinSystem.calculateTheoreticalMode(dir, qRange, structName, 'Type II', maxShear);
        end

        function [twin_objects, twin_structs] = calculateTheoreticalMode(inputObj, qRange, structName, type, maxShear)
            % Generic method to compute theoretical modes for both direct and reciprocal spaces.

            twin_objects_cell = {};          % Accumulate as cell to avoid subsref issues during sort
            twin_structs = struct.empty;     % Initialize empty struct array

            cs = inputObj.CS;
            isTypeI = strcmpi(type, 'Type I') || strcmpi(type, 'Type1') || strcmpi(type, 'type1');

            if isTypeI
                vals = [inputObj.h, inputObj.k, inputObj.l];
            else
                vals = [inputObj.u, inputObj.v, inputObj.w];
            end

            % Correct for floating point issues but throw error if any value is more than 1e-10 away from an integer
            if any(mod(vals, 1) > 1e-10 & mod(vals, 1) < (1 - 1e-10))
                error('Inputs must be integers (floating point noise excepted).');
            end

            v1 = round(vals(1)); v2 = round(vals(2)); v3 = round(vals(3));
            g = gcd(gcd(v1, v2), v3);

            % Enforce primitive plane/direction checks
            if g > 2
                twin_objects = twinSystem.empty;
                return;
            end

            % Apply selection rules
            if isTypeI
                if strcmpi(structName, 'FCC') && g == 1 && ~twinSystem.is_FCC_Rec(inputObj)
                    twin_objects = twinSystem.empty; return;
                end
                if strcmpi(structName, 'BCC') && g == 1 && ~twinSystem.is_BCC_Rec(inputObj)
                    twin_objects = twinSystem.empty; return;
                end
            else
                if strcmpi(structName, 'FCC') && g == 1 && ~twinSystem.is_FCC_Dir(inputObj)
                    twin_objects = twinSystem.empty; return;
                end
                if strcmpi(structName, 'BCC') && g == 1 && ~twinSystem.is_BCC_Dir(inputObj)
                    twin_objects = twinSystem.empty; return;
                end
            end
            if strcmpi(structName, 'none') && g > 1
                twin_objects = twinSystem.empty; return;
            end

            ap = [v1; v2; v3];

            if isTypeI
                [~, G_metric] = twinSystem.getMetricTensors(cs); % G_rec
            else
                [G_metric, ~] = twinSystem.getMetricTensors(cs); % G_dir
            end

            apD = G_metric * ap;
            d_sq = 1 / (ap' * apD);
            OH1_vec = d_sq * apD;

            % Solve the Diophantine equation to build valid basis vectors
            if isTypeI
                if g == 1
                    [OZ1, U, V] = twinSystem.Bezout3D(inputObj);
                    [~, U, V] = twinSystem.findUV(0.5*U, 0.5*V, structName);
                else
                    input_half = Miller(v1/2, v2/2, v3/2, cs, 'hkl');
                    [OZ1, U, V] = twinSystem.Bezout3D(input_half);
                    [oz2, U, V] = twinSystem.findUV(0.5*U, 0.5*V, structName);
                    OZ1 = Miller(0.5 * OZ1.u, 0.5 * OZ1.v, 0.5 * OZ1.w, cs, 'uvw');

                    if strcmpi(structName, 'FCC') && ~twinSystem.is_FCC_Dir(OZ1)
                        OZ1 = OZ1 + oz2;
                    elseif strcmpi(structName, 'BCC') && ~twinSystem.is_BCC_Dir(OZ1)
                        OZ1 = OZ1 + oz2;
                    end
                end
                BasisUV = [U.uvw', V.uvw'];
                OZ1_vec = OZ1.uvw';
            else
                if g == 1
                    [OZ1, U, V] = twinSystem.Bezout3DRec(inputObj);
                    [~, U, V] = twinSystem.findHK(0.5*U, 0.5*V, structName);
                else
                    input_half = Miller(v1/2, v2/2, v3/2, cs, 'uvw');
                    [OZ1, U, V] = twinSystem.Bezout3DRec(input_half);
                    [oz2, U, V] = twinSystem.findHK(0.5*U, 0.5*V, structName);
                    OZ1 = Miller(0.5 * OZ1.h, 0.5 * OZ1.k, 0.5 * OZ1.l, cs, 'hkl');

                    if strcmpi(structName, 'FCC') && ~twinSystem.is_FCC_Rec(OZ1)
                        OZ1 = OZ1 + oz2;
                    elseif strcmpi(structName, 'BCC') && ~twinSystem.is_BCC_Rec(OZ1)
                        OZ1 = OZ1 + oz2;
                    end
                end
                BasisUV = [U.hkl', V.hkl'];
                OZ1_vec = OZ1.hkl';
            end

            % Iterate through rational shear offsets
            for q = qRange

                % Stretch the single-layer theoretical normal (OH1_vec) by the twin index (q).
                % This point lies exactly vertically above the origin.
                OH = q * OH1_vec;

                % We stretch the single-layer physical atomic step (OZ1) by q.
                % This lands on an atom on the q-th plane, but because it stepped diagonally,
                % it is not necessarily the closest atom to the perfect normal (OH).
                OZ = q * OZ1_vec;

                % ZH points from our reference atom (OZ) to the geometric normal (OH).
                % Lies flat on the twin plane
                ZH = -OZ + OH;

                % This translates the 3D discrepancy (ZH) into our 2D (U, V) plane coordinate system.
                % How many fractional steps of U and V from reference atom to the geometric normal.
                ZH_inUV = BasisUV \ ZH;

                % Because atoms only exist at whole numbers, we find the bounding integers
                % of the fractional coordinates.
                nU = round(ZH_inUV(1));
                nV = round(ZH_inUV(2));

                rU = ZH_inUV(1) - nU;
                rV = ZH_inUV(2) - nV;

                sU = sign(rU); if sU == 0, sU = 1; end
                sV = sign(rV); if sV == 0, sV = 1; end

                % The 4 bounding nodes of the grid unit cell
                nodes = [
                    nU, nV;
                    nU + sU, nV;
                    nU, nV + sV;
                    nU + sU, nV + sV
                ];

                for nIdx = 1:4
                    cur_nU = nodes(nIdx, 1);
                    cur_nV = nodes(nIdx, 2);

                    % ZA is the 3D vector representing the flat physical translation along the twin plane
                    % from our reference atom (OZ) to this candidate real atom.
                    if isTypeI
                        ZA = cur_nU * U.uvw' + cur_nV * V.uvw';
                    else
                        ZA = cur_nU * U.hkl' + cur_nV * V.hkl';
                    end

                    % Add the adjustment walk (ZA) to our initial reference atom (OZ).
                    % OA is now the true vector from the origin to the candidate physical atom on the q-th plane.
                    % This is the shear direction (eta2).
                    OA = OZ + ZA;

                    % This is the remaining distance between our candidate physical atom (OA)
                    % and the geometric normal (OH).
                    AH = -ZA + ZH;

                    % OI is the vector from the origin to the sheared atom on the q-th plane.
                    ZI = ZA + 2 * AH;
                    OI = OZ + ZI;

                    % Validate elements and extract viable twin systems
                    eta2_vec = OA;
                    if norm(eta2_vec) < 1e-6
                        continue;
                    end

                    if isTypeI
                        vec2 = Miller(eta2_vec(1), eta2_vec(2), eta2_vec(3), cs, 'uvw');
                    else
                        vec2 = Miller(eta2_vec(1), eta2_vec(2), eta2_vec(3), cs, 'hkl');
                    end

                    % Build distortion matrix explicitly in lattice space
                    if isTypeI
                        BpOP = [U.uvw', V.uvw', OI];
                        BpOPd = [U.uvw', V.uvw', OA];
                    else
                        BpOP = [U.hkl', V.hkl', OI];
                        BpOPd = [U.hkl', V.hkl', OA];
                    end

                    if abs(det(BpOP)) < 1e-6
                        continue;
                    end

                    F_inner = BpOPd / BpOP;

                    if isTypeI
                        F = F_inner;
                    else
                        F = inv(F_inner)';
                    end

                    shearAmp = twinSystem.ShearMyFormula(F, cs);

                    % Filter and commit valid candidates
                    if shearAmp > 1e-4 && shearAmp <= maxShear
                        try
                            if isTypeI
                                tS = twinSystem.byK1eta2(inputObj, vec2);
                                rot = orientation.byAxisAngle(tS.k1, pi);
                            else
                                tS = twinSystem.byK2eta1(vec2, inputObj);
                                rot = orientation.byAxisAngle(tS.eta1, pi);
                            end

                            % Skip if the twin misorientation is a trivial symmetry operator
                            if angle(rot) < 1e-4
                                continue;
                            end

                            twin_objects_cell{end+1} = tS; %#ok<AGROW>

                            if isTypeI
                                twinStruct = struct('K1', inputObj, 'K2', tS.k2, 'eta1', tS.eta1, 'eta2', vec2, 'shear', shearAmp, 'q', q, 'F', F);
                            else
                                twinStruct = struct('K1', tS.k1, 'K2', vec2, 'eta1', inputObj, 'eta2', tS.eta2, 'shear', shearAmp, 'q', q, 'F', F);
                            end
                            twin_structs = [twin_structs; twinStruct]; %#ok<AGROW>
                        catch
                            % Skip if eta2 points directly along the plane nullifying shear
                            continue;
                        end
                    end
                end
            end

            % After iterating through all q values, sort the collected modes by shear magnitude
            if ~isempty(twin_structs)
                [~, sortIdx] = sort([twin_structs.shear]);
                twin_structs = twin_structs(sortIdx);

                twin_objects_cell = twin_objects_cell(sortIdx);
                twin_objects = cat(1, twin_objects_cell{:});
            else
                twin_objects = twinSystem.empty;
            end
        end
    end
end
