close all
clear all

logged = true;
fPlot = 408;

G0 = GlobalSkyModel2016;
G0 = G0.generate(fPlot);

%%
close all
G0.hideGround = false;

G1 = G0.underSample(3);
G1 = G1.setTime(datetime(2019,7,1,0,0,0));
% G1 = G1.setLocation([0,0,0]);
G1 = G1.setLocation('REACH');
C = CoordinateSystem;
% C = C.rotz(deg2rad(20));
% C = C.roty(deg2rad(20));
% C = C.rotx(deg2rad(20));
C = C.rotz(deg2rad(-90));
G1 = G1.setLocalCoor(C);
G1.projectionType = 'm';

figure,G1.plotProj(1,true)
G1.plotHorizon('k.')
G1.plotDirections('k')
G1.plotSun('k*')
G1.plotMoon('ko')
G1.plotVerifyMarkers('k')
G1.plotGridAz  
title('Galactic')

G2 = G1.changeGrid('RAdec');
figure,G2.plotProj(1,true)
G2.plotHorizon('k.')
G2.plotDirections('k')
G2.plotSun('k*')
G2.plotMoon('ko')
G2.plotVerifyMarkers('k')
G2.plotGridAz   
title('RADec')

G3 = G1.changeGrid('Horiz');
figure,G3.plotProj(1,true)
G3.plotHorizon('k.')
G3.plotDirections('k')
G3.plotSun('k*')
G3.plotMoon('ko')
G3.plotVerifyMarkers('k')
G3.plotGridAz
title('Horizontal')

G4 = G1.changeGrid('Local');
figure,G4.plotProj(1,true)
G4.plotHorizon('k.')
G4.plotDirections('k')
G4.plotSun('k*')
G4.plotMoon('ko')
G4.plotVerifyMarkers('k')
G4.plotGridAz
title('Local')

figure, G1.plotSkyView(1,true,[1,1,1,1])
title('Sky View')
figure, G4.plotSkyView(1,true,[1,1,1,1],true)
title('Sky View along local axis')

T3 = G3.interpOnHealPixGrid();
figure, healpixPlotMollweide(log2(T3))
title('Resampled on the Horizontal grid')
T4 = G4.interpOnHealPixGrid();
figure, healpixPlotMollweide(log2(T4))
title('Resampled on the Local grid')

figure
G4.plotLocalCoor


