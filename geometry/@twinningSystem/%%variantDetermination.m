%function [isTwinning, variantIndex,k1, rotAxis] = variantDetermination(tS, misOrientation, deviationInDegree)
%            
%            
%tSvariants = tS.symmetrise('antipodal'); 
%            
%variantMisorientations = arrayfun(@(x) orientation.map(x.n,x.n,x.b,-x.b), tSvariants);
%
%            
%variantMatch = angle(misOrientation, transpose(variantMisorientations),'noSym1')/degree;
%            
%create logical array for match between variant misorientation and grain misorientation
%            
%variantMatch_logic= variantMatch==min(variantMatch,[],2)&variantMatch<deviationInDegree;
%            
%isTwinning = logical(variantMatch_logic*ones(length(tSvariants),1));
%            
%variantIndex = variantMatch_logic*transpose(1:length(tSvariants));
%            
%variantIndex = variantIndex(isTwinning);%
%
%tSvariants(2)
%
%vars = tSvariants(variantIndex);
%
%k1 = vars.n;
%
%rotAxis = vars.b;
%
%end