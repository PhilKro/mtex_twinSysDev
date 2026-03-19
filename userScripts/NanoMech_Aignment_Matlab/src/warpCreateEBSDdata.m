function [newEBSD] = warpCreateEBSDdata(oldEBSD,Xnew,Ynew,alignmentData)
%%% this function creates a new Mtex EBSD variable/struct by warping an old one using alignment data created from point selection.
%%% Xnew and Ynew are the aimed grids in the target configuration (i.e.
%%% reference), align_data contains the warping information

%%% author: Tijmen Vermeij, TUe, 2018-2023

if any(size(oldEBSD.x)==1)
    error('ebsd input data must be gridded')
end


if nargin == 4
    % get relevant displacement fields by interpolating from thez original to the new grid of positions in
    % the ref config
    Ux = interp2(alignmentData.X,alignmentData.Y,alignmentData.Ux,Xnew,Ynew);
    Uy = interp2(alignmentData.X,alignmentData.Y,alignmentData.Uy,Xnew,Ynew);
else
    warning('no alignmentData given, so no warping is performed, but only interpolation.')
    Ux = zeros(size(Xnew));
    Uy = zeros(size(Xnew));
end

% use MTEX interp function to interpolate the old EBSD data to the new
% [newEBSD, idNearest] = interp(oldEBSD,Xnew+Ux,Ynew+Uy,'method','nearest');
newEBSD = interp(oldEBSD,vector3d(Xnew+Ux,Ynew+Uy,0),'method','nearest');
newEBSD.x = Xnew(:);
newEBSD.y = Ynew(:);

% create new unitCell variable, to be sure.
unitCellXY = calcUnitCell([Xnew(:) - mean(Xnew(:)),Ynew(:) - mean(Ynew(:))]);
warning('warpCreateEBSDdata: this assignment will break at some point when calcUnitCell is fixed to returen a vector3d')
newEBSD.unitCell = vector3d(unitCellXY(:,1),unitCellXY(:,2), zeros(size(unitCellXY(:,2))));

% idNearest is not returned by interp, at least not in MTEX 6.1.1
% % return idNearest as varargout
% if nargout > 1
%     varargout{1} = idNearest;     % <<< added
% end

end