function [tSsym,id] = symmetrise(tS,varargin)
% find all symmetrically equivalent twin systems
%
% Syntax
%
%   tSAll = tS.symmetrise
%   [tSAll,id] = symmetrise(tS)
%
% Input
%  tS - @twinSystem
%
% Output
%  tS - @twinSystem
%  id    - id of the twinSystem before symmetrisation
%
%% 


% find all symmetrically equivalent twin systems

if ~isa(tS.eta1,'Miller'), id = 1:length(tS); return; end

id = [];

for i = 1:length(tS)
    
    s = tS(i).CS.quaternion; % all proper symmetry operations
    
    symcheck = s * tS(i).eta1;

    if ~(check_option(varargin,'antipodal')) 
        [~, unique_idx] = unique(symcheck, 'antipodal', 'noSymmetry', 'stable');
    else
        % If the flag is true, keep all indices
        % We use a column vector to match the typical output shape of 'unique'
        unique_idx = (1:length(symcheck))'; 
    end
    
    
    % apply symmetry operations to all twinning vectors
    eta1_s = s(unique_idx) * tS(i).eta1;
    k1_s = s(unique_idx) * tS(i).k1;
    k2_s = s(unique_idx) * tS(i).k2;
    
    % Sanity check that all planes still have ths same relation with each
    % other
    tol = 1e-4; 
    ang_k1_eta1 = angle(k1_s, eta1_s, 'noSymmetry');
    ang_k2_eta1 = angle(k2_s, eta1_s, 'noSymmetry');
    ang_k1_k2   = angle(k1_s, k2_s, 'noSymmetry');
    check1 = (max(ang_k1_eta1) - min(ang_k1_eta1)) < tol;
    check2 = (max(ang_k2_eta1) - min(ang_k2_eta1)) < tol;
    check3 = (max(ang_k1_k2)   - min(ang_k1_k2))   < tol;
    if ~(check1 && check2 && check3)
        error('K2, K1 and eta1 dont have the same relation with each other in the symmetrised objects.');
    end
    
    % now we have the unique symmetrised triplets, we can create the twinSystem objects
    
    rotAxis_s = cross(eta1_s,k1_s);
    eta2_s = cross(rotAxis_s,k2_s);

    if tS(i).variantId == 1
        variantId_s = transpose(1:length(eta1_s));
    else
        warning('twinSystem:symmetrise:notVariant1','Hope the variant assignment works correctly. Please check the output carefully.')
        variantId_s = transpose(mod((tS(i).variantId - 1) + (0:length(eta1_s)-1), length(eta1_s)) + 1);
    end

    
    tS_new = twinSystem(rotAxis_s, k1_s, eta1_s, k2_s, eta2_s, repmat(tS(i).CRSS,length(eta1_s),1), repmat(tS(i).twinType,length(eta1_s),1), variantId_s);
    
    if i == 1
        tSsym = tS_new;
    else
        tSsym = [tSsym; tS_new];
    end
    
    id = [id; repmat(i, length(tS_new), 1)];
end

end


% if ~isa(tS.eta1,'Miller'), return; end

% CRSS = [];
% id = [];
% k1 = [];
% k2 = [];
% eta1 = [];

% for i = 1:length(tS)

%     % find all symmetrically equivalent
%     if tS.CS.lattice == 6
%         %easier like this?
%         %orientation.byAxisAngle(Miller(0,0,0,1,CS),60*degree, CS,CS)
%         % antipodal 
%         % -orientation.byAxisAngle(Miller(0,0,0,1,CS),180*degree, CS,CS)




%         % a "directional antipodal" keeps the same sense with respect to the
%         % c-axis
%         % In particular, twinning induces a well-defined constant amount of
%         % shear and is intrinsically polarized. The latter means that, 
%         % unlike the plastic slip, the activation of twin depends on the 
%         % sense of the shear within the twin plane.
%         % (https://doi.org/10.1016/j.jmps.2022.104855)
%         mm = symmetrise(tS.eta1(i),'unique');
%         mm = mm(isnull(dot(mm, Miller({0,0,0,1},tS.CS),'noSymmetry')-dot(tS.eta1(i), Miller({0,0,0,1},tS.CS),'noSymmetry')));
        
%         nn = symmetrise(tS.k1(i),'unique'); %#ok<*PROP>
%         nn = nn(isnull(dot(nn, Miller({0,0,0,1},tS.CS),'noSymmetry')-dot(tS.k1(i), Miller({0,0,0,1},tS.CS),'noSymmetry')));
%         % find those which have the same angles as the original system
%         % for twinning System K1 (n) and the eta1 (b) that is 90°
%         [r,c] = find(isnull(dot(tS.k1(i),tS.eta1(i),'noSymmetry')-dot_outer(mm,nn,'noSymmetry')));

%         % restricht to the orthogonal ones
%         eta1 = [eta1;mm(r(:))]; %#ok<*AGROW>
%         k1 = [k1;nn(c(:))];

%         % find the corresponding k2 planes under the same angle to k1
%         oo = symmetrise(tS.k2(i),'unique');
%         oo = oo(abs(dot(oo, Miller({0,0,0,1},tS.CS),'noSymmetry')-dot(tS.k2(i), Miller({0,0,0,1},tS.CS),'noSymmetry'))<1.00e-012);

%         % find k2s which have the same angles between k2 and K1 (n)
%         [~, rk2] = find(isnull(dot(tS.k2(i),tS.k1(i),'noSymmetry')-dot_outer(nn,oo,'noSymmetry')));
%         %[rk2, ~] = find(isnull(dot(tS.n(i),tS.k2(i),'noSymmetry')-dot_outer(n,k2s,'noSymmetry')));
%         %PROBLEM BASAL PLANE AS K2
%         k2 = [k2; oo(rk2(:))];

%         rotAxis = cross(eta1,k1);
%         eta2 = cross(rotAxis,k2);
%         rotAxis.dispStyle='UVTW';

%         %Checks
%         assert(all(isnull(dot(k1,k2,'noSymmetry')-dot(tS.k1(i),tS.k2(i),'noSymmetry'))),'Angles between K1 and K2 are not the same')
%         assert(all(isnull(dot(k1,eta1,'noSymmetry'))), 'Not all eta1s lie in k1s')
%         shearvec = 2*(k1.normalize-dot(eta2.normalize,k1.normalize).^-1.*eta2.normalize);
%         assert(all(dot(shearvec, eta1,'noSymmetry'))>0, "Something went wrong in the definition of the twin system: eta1 has the wrong sense. Sense of shear will not be correct")
%     else
%         %incomplete
%         mm = symmetrise(tS.eta1(i),'unique','antipodal');
%         if ~check_option(varargin,'antipodal'), mm = [mm;-mm]; end
%         nn = symmetrise(tS.k1(i),'unique','antipodal'); %#ok<*PROP>
%         oo = symmetrise(tS.k2(i),'unique','antipodal');
%         %DOES THIS MAKE SENSE
%         if ~check_option(varargin,'antipodal'), oo = [oo;-oo]; end
%     end


%     CRSS = [CRSS;repmat(tS.CRSS(i),length(r),1)];
%     id = [id;repmat(i,length(r),1)];

%     % not yet defined
% end



% tS.k1 = k1;
% tS.rotAxis = rotAxis;
% tS.eta1 = eta1;
% tS.k2 = k2;
% tS.eta2 = eta2;

% tS.CRSS = CRSS;
% end
