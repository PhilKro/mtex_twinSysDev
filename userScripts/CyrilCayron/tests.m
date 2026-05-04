CS = crystalSymmetry('6/mmm',[3.21 3.21 5.21],'x||a','mineral','Magnesium');

% [testTwins, testStruct] =  twinSystem.calculateTheoreticalTwins(CS,2,1);
% [testTwins, testStruct] =  twinSystem.calculateTheoreticalTwins(CS,Miller(0,1,-1,2,CS),8);
[testTwins, testStruct] =  twinSystem.calculateTheoreticalTwins(CS,Miller(1,1,-2,2,CS),8);
% [testTwins, testStruct] =  twinSystem.calculateTheoreticalTwins(CS,Miller(1,1,-2,1,CS),8);
% [testTwins, testStruct] =  twinSystem.calculateTheoreticalTwins(CS,Miller(0,0,0,1,CS),8);
twinSystem.printTheoreticalTwins(testStruct);
testTwins

%% cubic FCC
CS_cube = crystalSymmetry('m3m',[1 1 1],'mineral','cubic');
[testTwins, testStruct] =  twinSystem.calculateTheoreticalTwins(CS_cube,Miller(1,1,1,CS_cube),8, 'FCC');

twinSystem.printTheoreticalTwins(testStruct);
testTwins
%% cubic BCC
CS_cube = crystalSymmetry('m3m',[1 1 1],'mineral','cubic');
[testTwins, testStruct] =  twinSystem.calculateTheoreticalTwins(CS_cube,Miller(1,1,2,CS_cube),8, 'BCC');

twinSystem.printTheoreticalTwins(testStruct);
testTwins
%% plot shear plane

twinSystem.plotShearPlane(Miller(0,1,-1,2,CS), 1:5)
% SumTingWong