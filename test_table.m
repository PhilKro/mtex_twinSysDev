T = table([1; 2], {'a'; 'b'}, 'VariableNames', {'Id', 'Val'});
try
    val = T.Val{0};
    disp('No crash: '); disp(val);
catch ME
    disp('Crashed: '); disp(ME.message);
end
quit;
