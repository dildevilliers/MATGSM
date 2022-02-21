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
G1 = G1.setTime(datetime(2021,09,20,12,0,0));
G1 = G1.setLocation([0,0,0]);
G1.projectionType = 'm';

figure,G1.plotProj(1,true)
G1.plotHorizon('k.')
G1.plotDirections('k')
G1.plotSun('k*')
G1.plotMoon('ko')
G1.plotVerifyMarkers('k')

G2 = G1.changeGrid('RAdec');
figure,G2.plotProj(1,true)
G2.plotHorizon('k.')
G2.plotDirections('k')
G2.plotSun('k*')
G2.plotMoon('ko')
G2.plotVerifyMarkers('k')

G3 = G1.changeGrid('Horiz');
figure,G3.plotProj(1,true)
G3.plotHorizon('k.')
G3.plotDirections('k')
G3.plotSun('k*')
G3.plotMoon('ko')
G3.plotVerifyMarkers('k')

figure, G1.plotSkyView(1,true,[1,1,1,1])

