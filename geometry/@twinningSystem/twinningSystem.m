classdef twinningSystem < slipSystem
    properties
        eta1
        eta2
        k1
        k2
        rotAxis % better name: zoneAxis
    end

    methods
        function tS = twinningSystem(rotAxis, k1, eta1, k2, eta2)
            tS = tS@slipSystem(rotAxis, k1);
            tS.k1 = tS.n;
            
            % THIS HAS TO BE CHANGED TO b = eta 1 . repurcussions on the
            % symmetrise functions
            tS.eta1 = eta1;
            tS.rotAxis = tS.b;
            tS.k2 = k2;
            tS.eta2 = eta2;
        end

        function [isTwinning, variantIndex,tSvariants] = variantDetermination(tS, misOrientation, deviationInDegree)
            tSvariants = tS.symmetrise('antipodal');
            variantMatch = angle(misOrientation, tSvariants.parentTwinMisorientation(),'noSym1')/degree;
            %variantMatch = cell2mat(arrayfun(@(x) angle(misOrientation, x,'noSym1')/degree, tSvariants.parentTwinMisorientation(), 'UniformOutput', false));
            %create logical array for match between variant misorientation and grain misorientation
            variantMatch_logic= variantMatch==min(variantMatch,[],2)&variantMatch<deviationInDegree;
            isTwinning = logical(variantMatch_logic*ones(length(tSvariants),1));
            variantIndex = variantMatch_logic*transpose(1:length(tSvariants));
            variantIndex = variantIndex(isTwinning);
        end
        
        function [variantIndex] = variantDetermination2(tS, misOrientation, deviationInDegree)
            % Input
            %  misOrientation - n x 1 array of misOrientations
            %  deviationInDegree - double, allowed error for twin
            %    determination
            %
            % Output
            %  variantIndex - n x size(tS) double, where 0 = no twin system; 
            %  m = (1 : length(tS.symmetrise)) = indicates the variant

            tSvariants = tS.symmetrise('antipodal');
            %angle with 'noSym1' needs both column or both row vectors
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
            variantMatch_angle = angle(misOrientation, variantPTM,'noSym1')/degree;
            %variantMatch_angle = cell2mat(arrayfun(@(x) angle(misOrientation, x,'noSym1'), tSvariants.parentTwinMisorientation(), 'UniformOutput', false))./degree;
            %create logical array for match between variant misorientation and grain misorientation
            variantMatch_logic= variantMatch_angle==min(variantMatch_angle,[],2)&variantMatch_angle<deviationInDegree;
            %isTwinning = variantMatch_logic*ones(length(tSvariants),1)
            %%isTwinning = logical(variantMatch_logic*ones(length(tSvariants),1));
            
            % this works if misori is defined by: misori = orientation.map(tS.n,tS.n,tS.b,-tS.b);
            variantIndex = variantMatch_logic*transpose(1:length(tSvariants));
            anglesMisoriFromTwin = min(variantMatch_angle, [], 2);

            % this works if misori is defined by improper rotation: misori = -orientation.byAxisAngle(tS.k1, pi);
            % variantIndex = transpose((1:length(tSvariants))*variantMatch_logic);
        end
        


        function [anglesMisoriFromTwin] = anglefromTwin(tS, misOrientation)
            %this is bullshit (should be a second output from variantDetermination2) but i dont get two outputs through the subsref method
            % Input
            %  misOrientation - n x 1 array of misOrientations
            %  deviationInDegree - double, allowed error for twin
            %    determination
            %
            % Output
            %  variantIndex - n x size(tS) double, where 0 = no twin system; 
            %  m = (1 : length(tS.symmetrise)) = indicates the variant

            tSvariants = tS.symmetrise('antipodal');
            %angle with 'noSym1' needs both column or both row vectors
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
            variantMatch_angle = angle(misOrientation, variantPTM,'noSym1')/degree;
            anglesMisoriFromTwin = min(variantMatch_angle, [], 2);
        end

        function misori = parentTwinMisorientation(tS)
            % parentTwinMisorientation
            %
            % Syntax
            %  misori = parentTwinMisorientation(tS)
            %
            % Input
            %  tS - list of @twinningSystem
            %
            % Output
            %  misori - misorientation
            %  in crystal Coordinates parentOrientataion *
            %  parentTwinMisorientation == parentOrientataion *
            %  inv(parentTwinMisorientation)
            %  bc its a mirroriong (only the quaterinons flip the sign)
            misori = orientation.map(tS.n,tS.n,tS.b,-tS.b);
            
            % https://mtex-toolbox.github.io/RotationImproper.html
            % results in the same euler angles jus the Inv. = true
            % misori.i = true;
            % misori = -orientation.byAxisAngle(tS.k1, pi);
            assert(max(angle(misori,-orientation.byAxisAngle(tS.k1, pi)),[],'all')<1e-5, "Mapping of n and b is not equivalent to mirroring across K1")
        end

        function shear = shearMagnitude(tS)
            shear = sqrt(4*((dot(tS.eta2.normalize,tS.k1.normalize).^-2)-1));
        end
        
        function [shearvector, shearMagnitude] = shear(tS)
            shearvector = 2*(tS.k1.normalize-dot(tS.eta2.normalize,tS.k1.normalize).^-1*tS.eta2.normalize);
            shearvector.dispStyle='UVTW';
            shearMagnitude = norm(shearvector);
            %shear = sqrt(4*((dot(tS.eta2.normalize,tS.k1.normalize).^-2)-1));
        end
        
        function defTensor = displacementGradient(tS)
            % returns deformation tensor in parent reference frame
            % Deformation Gradient tensor in specimen geometry 
            % = tensor.eye + parentOrientation * tS.deformationTensor
            defTensor = (tS.shearMagnitude.*dyad(tS.eta1, tS.k1).normalize);
        end
        
        function correspondanceMatrix = correspondanceMatrix(tS)
            correspondanceMatrix = -tensor.eye(tS.CS) + 2*dyad(tS.eta2, tS.k1).normalize;
        end

        function nrVariants = numberOfVariants(tS)
            nrVariants = length(tS.symmetrise('antipodal'));
        end

        function display(tS,varargin) %#ok<DISPLAY>
            % standard output

            displayClass(tS,inputname(1),varargin{:},'moreInfo',char(tS.CS,'compact'));

            if length(tS)>24, disp([' CRSS: ' xnum2str(unique(tS.CRSS))]); end
            if length(tS)>1, disp([' size: ' size2str(tS.b)]); end

            disp(' ');

            if length(tS)<=45 && ~isempty(tS)
                dispData(tS)
            elseif ~getMTEXpref('generatingHelpMode')
                disp(' ')
                setappdata(0,'dispLastData',@() dispData(tS));
                disp('  <a href="matlab:feval(getappdata(0,''dispLastData''))">display all coordinates</a>')
                disp(' ')
            end
        end
        function n = char(tS)
            reta1 = char(round(tS.eta1));
            rk1 = char(round(tS.k1));
            n = strcat(rk1, " ",reta1);
        end
        function n = filename_char(tS)
            reta1 = char(round(tS.eta1));
            rk1 = char(round(tS.k1));
            n = strcat("p",rk1(2:end-1), "d",reta1(2:end-1));
            n = strrep(n, '-', '_');
        end
        function dispData(tS)
            % display coordinates
            if isa(tS.CS,'crystalSymmetry')
                if tS.eta1.lattice.isTriHex
                    reta1 = round(tS.eta1);
                    rk1 = round(tS.k1);
                    d = [reta1.UVTW rk1.hkil];
                    d(abs(d) < 1e-10) = 0;
                    cprintf([d,reshape(tS.CRSS,[],1)],'-L','  ','-Lc',{'eta_1 U' 'V' 'T' 'W' '| K_1 H' 'K' 'I' 'L' 'CRSS'});
                else
                    d = [tS.eta1.uvw tS.k1.hkl];
                    d(abs(d) < 1e-10) = 0;
                    cprintf([d,reshape(tS.CRSS,[],1)],'-L','  ','-Lc',{'eta_1 u' 'v' 'w' '| K_1 h' 'k' 'l' 'CRSS'});
                end
            else
                d = round(100*[tS.eta1.xyz tS.k1.xyz])./100;
                d(abs(d) < 1e-10) = 0;
                cprintf(d,'-L','  ','-Lc',{'x' 'y' 'z' ' |   x' 'y' 'z' });
            end
        end
    end

    methods (Static = true)
        function tS = byK1eta2(k1, eta2)
            rotAxis = cross(k1, eta2);
            eta1 = cross(k1, rotAxis);
            k2 = cross(eta2, rotAxis);
            tS = twinningSystem(rotAxis, k1, eta1, k2, eta2);
            %checks
            shearvec = 2*(tS.k1.normalize-dot(tS.eta2.normalize,tS.k1.normalize).^-1*tS.eta2.normalize);
            assert(dot(shearvec, tS.eta1,'noSymmetry')>0, "Something went wrong in the definition of the twin system: eta1 has the wrong sense. Sense of shear will not be correct")
        end

        function tS = byK2eta1(k2, eta1)
            rotAxis = cross(eta1,k2);
            rotAxis.dispStyle='UVTW';
            k1 = cross(rotAxis, eta1);
            k1.dispStyle = 'hkil';
            eta2 = cross(rotAxis,k2);
            eta2.dispStyle='UVTW';
            tS = twinningSystem(rotAxis, k1, eta1, k2, eta2);
            %checks
            shearvec = 2*(tS.k1.normalize-dot(tS.eta2.normalize,tS.k1.normalize).^-1*tS.eta2.normalize);
            assert(dot(shearvec, tS.eta1,'noSymmetry')>0, "Something went wrong in the definition of the twin system: eta1 has the wrong sense. Sense of shear will not be correct")
        end

        function tS = byK1rotAxis(k1, rotAxis)
            if dot(k1, rotAxis,'noSymmetry')<1.00e-012
                eta1 = cross(k1,rotAxis);
                eta1.dispStyle='UVTW';
                k2 = eta1;
                k2.dispStyle = 'hkil';
                eta2 = k1;
                eta2.dispStyle='UVTW';
                disp(append(newline,'Chose a set of values for K2 and eta2!'))
                tS = twinningSystem(rotAxis, k1, eta1, k2, eta2);
            else
                disp("k1 and rotAxis have to be orthogonal")
            end
        end

        function tS = hexagonal_110K(K, CS)
            %constructs twinningSystem objects for {1, 0 -1, K} type twins
            %most commonly K ∊ [1 2 3]
            %(https://journals.aps.org/prb/abstract/10.1103/PhysRevB.78.024117)
            % Input
            % K - array of float
            % CS - hexagonal crystal symmetry
            if CS.lattice == 6
                tS = arrayfun(@(x_int) twinningSystem.byK1rotAxis(Miller({1,0,-1,x_int}, CS), Miller({1,-2,1,0}, CS, 'UVTW')), K);
                % twinningSystem.byK1rotAxis(Miller({1,0,-1,K}, CS), Miller({1,-2,1,0}, CS, 'UVTW'));
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_211K(K, CS)
            %constructs twinningSystem objects for {1, 1, -2, K} type twins
            %most commonly K ∊ [1 2 3 4]
            %(https://journals.aps.org/prb/abstract/10.1103/PhysRevB.78.024117)
            % Input
            % K - array of float
            % CS - hexagonal crystal symmetry

            if CS.lattice == 6
                tS = arrayfun(@(x_int) twinningSystem.byK1rotAxis(Miller({1,1,-2,x_int}, CS), Miller({1,-1,0,0}, CS, 'UVTW')), K);
                %tS = twinningSystem.byK1rotAxis(Miller({1,1,-2,K}, CS), Miller({1,-1,0,0}, CS, 'UVTW'));
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1012_1012(CS)
            if CS.lattice == 6
                tS = twinningSystem.byK2eta1(Miller(-1,0,1,2,CS), Miller(-1,0,1,1,CS,'UVTW'));
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1122_0001(CS)
            %opposite eta1 direction to hexagonal_1122_1124
            if CS.lattice == 6
                tS = twinningSystem.byK2eta1(Miller(0,0,0,1,CS), Miller(-1,-1,2,3,CS,'UVTW'));
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1101_1013(CS)
            %opposite eta1 direction to hexagonal_1122_1124
            if CS.lattice == 6
                tS = twinningSystem.byK2eta1(Miller(1,0,-1,-3,CS), Miller(1,0,-1,-2,CS,'UVTW'));
            else
                disp("Provide a hexagonal CS for this method")
            end
        end

        function tS = hexagonal_1122_1124(CS)
            %opposite eta1 direction to hexagonal_1122_0001
            if CS.lattice == 6
                tS = twinningSystem.byK2eta1(Miller(1,1,-2,-4,CS), Miller(1,1,-2,-3,CS,'UVTW'));
            else
                disp("Provide a hexagonal CS for this method")
            end
        end
        
        function tS = hexagonal_1124_1122(CS)
            if CS.lattice == 6
                tS = twinningSystem.byK2eta1(Miller(1,1,-2,-2,CS), Miller(2,2,-4,-3,CS,'UVTW'));
            else
                disp("Provide a hexagonal CS for this method")
            end
        end


        function tS = hexagonal_1121_0001(CS)
            if CS.lattice == 6
                tS = twinningSystem.byK2eta1(Miller(0,0,0,1,CS), Miller(-1,-1,2,6,CS,'UVTW'));
            else
                disp("Provide a hexagonal CS for this method")
            end
        end
    end
end

