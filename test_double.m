a = 1;
a(1) = [];
disp(size(a));
try
   a = [a; 2];
   disp('Success double vertcat');
catch ME
   disp(ME.message);
end

try
   a = [a, 2];
   disp('Success double horzcat');
catch ME
   disp(ME.message);
end
quit;
