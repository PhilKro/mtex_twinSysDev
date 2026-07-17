ori = orientation.id(crystalSymmetry('m-3m'));
sq_centers = ori;
sq_centers(1) = [];
disp('Size after deletion:');
disp(size(sq_centers));
try
   sq_centers = [sq_centers; ori];
   disp('Success vertcat');
catch ME
   disp('Error vertcat:');
   disp(ME.message);
end

try
   sq_centers = [sq_centers, ori];
   disp('Success horzcat');
catch ME
   disp('Error horzcat:');
   disp(ME.message);
end

% We can also just reshape or create an empty correctly
sq_empty = orientation.empty(0, 1);
try
   sq_empty = [sq_empty; ori];
   disp('Success vertcat with empty(0,1)');
catch ME
   disp('Error vertcat empty(0,1):');
   disp(ME.message);
end
quit;
