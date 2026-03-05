function c = mtimes(a,b)
% implement multiplication for grain2d

if isa(b,'twinSystem')
  c = activate(b,a);
  return
end

error('mtex:grain2d:mtimes','Undefined function mtimes for grain2d.');
end
