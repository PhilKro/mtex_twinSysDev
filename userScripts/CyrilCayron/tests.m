CS = crystalSymmetry('6/mmm',[3.21 3.21 5.21],'x||a','mineral','Magnesium');

[testTwins, testStruct] =  twinSystem.calculateTheoreticalTwins(CS,2,6);
angle(testTwins.k1, Miller(1,0,-1,2,testTwins.CS))<1e-5;
log1012 = angle(testTwins.k1, Miller(0,1,-1,2,testTwins.CS),'noSymmetry')<1e-5;

twinSystem.printTheoreticalTwins(testStruct(log1012))
% COMPARE TO THE ORIGINAL PROGRAM QUITE SOME DIFFERENCES, SOME STUFF IS CORRECT SOME NOT