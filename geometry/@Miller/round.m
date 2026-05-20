function h = round(h,varargin)
% tries to round miller indizes to greatest common divisor
%
% Syntax
%   h = round(h)
%   h = round(h,'maxHKL',5)
%   h = round(h,'use4Index')
%
% Input
%  h - @Miller
%
% Output
%  h - @Miller
%
% Options
%  maxHKL - maximum value of Miller indices
%  use4Index - flag to decide if all 4 indices should be used for rounding
%

% ignore xyz case
if h.dispStyle == MillerConvention.xyz, return; end

sh = size(h);

mOld = h.coordinates;

% check if 4 index notation should be used
use4Index = check_option(varargin, 'use4Index');

% consider only 3 digits Miller indices unless use4Index is specified
if use4Index && size(mOld, 2) == 4
  mOld = mOld';
else
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

% remember current dispStyle
currentDispStyle = h.dispStyle;

% now round
h.coordinates = round(h.coordinates);

% restore dispStyle to prevent breaking it (e.g. hkl changing to hkil)
h.dispStyle = currentDispStyle;

h = reshape(h,sh);
