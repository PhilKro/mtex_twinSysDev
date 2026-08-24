CS = crystalSymmetry('6/mmm');
ori = orientation.rand(7, CS);
colors = rand(7, 3);
figure;
try
    plotPDF(ori, Miller(0,0,0,1,CS), 'property', colors);
    disp('Syntax 1 works: property');
catch ME
    disp(['Syntax 1 failed: ', ME.message]);
end

figure;
try
    plotPDF(ori, colors, Miller(0,0,0,1,CS));
    disp('Syntax 2 works: positional colors');
catch ME
    disp(['Syntax 2 failed: ', ME.message]);
end

figure;
try
    % Forcing scatter by using the scatter function? plotPDF handles symmetries.
    plot(ori * Miller(0,0,0,1,CS), 'MarkerFaceColor', colors);
    disp('Syntax 3 works: plot with projected vectors');
catch ME
    disp(['Syntax 3 failed: ', ME.message]);
end
