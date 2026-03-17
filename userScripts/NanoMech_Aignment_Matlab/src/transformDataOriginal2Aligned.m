function aligned_field = transformDataOriginal2Aligned(field, originalEBSD, ebsd)
%ALIGNGND_E Maps the field data vector based on the 'originalId' property
% from 'originalEBSD' dataset to the aligned 'ebsd' dataset.
%
% Syntax:
%   aligned_GND_E = alignGND_E(GND_E, originalEBSD, ebsd)
%
% Inputs:
%   field          - A column vector of double data (e.g., geometrically necessary 
%                    dislocation densities) corresponding to the data points 
%                    in originalEBSD.
%   originalEBSD   - The original EBSD dataset (MTEX-style variable) 
%                    containing the 'originalId' property.
%   ebsd           - The modified or subset EBSD dataset (MTEX-style 
%                    variable) for which the data needs to be aligned.
%
% Output:
%   aligned_GND_E  - A column vector of the same length as 'ebsd', 
%                    where each entry corresponds to the GND_E value of the
%                    data point in 'originalEBSD' that matches the 
%                    'originalId' of the corresponding point in 'ebsd'.
%
% Dependencies: This function requires MTEX and assumes the MTEX variables
%               have a 'prop' struct with 'originalId' fields and the 
%               prop2ind method (modified version of id2ind).

    % Get the size of the input data vector (GND_E)
    sW = size(field);

    % 1. Create a version of GND_E with an extra NaN row at the end.
    % This NaN entry will be used for EBSD points that don't map back
    % to the original data (where listInds will be 0).
    % Ensure GND_E is a column vector before concatenation if it isn't
    if sW(2) > sW(1)
        field = field';
        sW = size(field);
    end
    
    GND_E_withNanRow = [field; NaN(1, sW(2))];
    lastNanEntryId = length(GND_E_withNanRow);

    % 2. Extract the indices for a given originalId which point to the same data
    % in originalEBSD (i.e., originalEBSD(listInds) corresponds to ebsd).
    % The ebsd.prop.originalId is used as the key for the mapping.
    listInds = originalEBSD.prop2ind('originalId', ebsd.prop.originalId);

    % 3. Where listInds are 0 (meaning no match was found in originalEBSD
    % for an 'ebsd' point's 'originalId'), default the index to the 
    % ID of the last entry, which is the NaN value.
    listInds(listInds == 0) = lastNanEntryId;

    % 4. Use the aligned indices to create the new vector.
    aligned_field = GND_E_withNanRow(listInds);

end