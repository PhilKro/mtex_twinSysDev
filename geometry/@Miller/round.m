function [h, err] = round(h,varargin)
% tries to round miller indizes to greatest common divisor
%
% Syntax
%   h = round(h)
%   h = round(h,'maxHKL',5)
%   h = round(h,'round3Index')
%   [h, err] = round(...)
%
% Input
%  h - @Miller
%
% Output
%  h - @Miller
%  err - angle between original and rounded Miller vector
%
% Options
%  maxHKL - maximum value of Miller indices
%  round3Index - flag to force rounding the 3 index notation for Tri/Hex lattices
%

% ignore xyz case
if h.dispStyle == MillerConvention.xyz, err = 0; return; end

h_old = h;
sh = size(h);

% check if we should force 3 index rounding for Tri/Hex
round3Index = check_option(varargin, 'round3Index') && h.lattice.isTriHex;

isDirection = any(strcmpi(char(h(1).dispStyle), {'uvw','UVTW'}));

if round3Index
  if isDirection
      mOld = h.uvw;
  else
      mOld = h.hkl;
  end
  mOld = mOld';
else
  mOld = h.coordinates;
  % consider only 3 digits Miller indices
  mOld = mOld(:,[1 2 end])';
end

% the 
mMax = reshape(max(abs(mOld),[],1),size(h));
%mbr = reshape(selectMaxbyColumn(abs(mv)),size(h));

maxHKL = get_option(varargin,'maxHKL',12);

multiplier = ones(size(h));
for im = 1:size(mOld,2)
  
  mNew = mOld(:,im) / mMax(im) * (1:maxHKL);  
  
  e = 1e-7 * round(1e7 * sum((mNew - round(mNew)).^2)./sum(mNew.^2));
    
  [~,n] = min(e);
  
  multiplier(im) = n / mMax(im);
  
end

h = h .* multiplier;

currentDispStyle = h.dispStyle;

% now round
if round3Index
  if isDirection
      h.uvw = round(h.uvw);
  else
      h.hkl = round(h.hkl);
  end
else
  h.coordinates = round(h.coordinates);
end

h.dispStyle = currentDispStyle;

h = reshape(h,sh);

if nargout == 2
    err = angle(h_old, h);
end
