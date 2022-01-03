close all
clear all

logged = true;

% GSM = GlobalSkyModel([],'haslam');
% GSM = GSM.generate([408]);
% GSM.view(1,logged)
% GSM.write_fits('c:\Temp\GSMfits.fits')
% 
% %%
% signPhi = -1;
% map = GSM.generated_map_data(:,1);
% if logged, map = log2(map); end
% sz = size(map);
% nside = sqrt(max(sz)/12);
% tp = pix2ang(nside);
% tp = [tp{:}];
% th = tp(1,:);
% ph = signPhi.*tp(2,:);
% lat = pi/2 - th(:);
% long = wrap2pi(ph(:));
% 
% figure
% gridDelaunay = delaunay(long,lat);
% h = trisurf(gridDelaunay,rad2deg(long),rad2deg(lat),map.*0,map);
% set(h, 'LineStyle', 'none')
% axis equal
% % axis off
% view([0,90])
% xlabel('Long^\circ')
% ylabel('Lat^\circ')
% 
% 
% time = datetime(2018,7,22,0,0,0);
% julDate = convert.date2jd([time.Day,time.Month,time.Year,time.Hour,time.Minute,time.Second]);
% % earthLocation = [deg2rad(18.86) deg2rad(-33.93) 300];
% earthLocation = [deg2rad(25) deg2rad(0) 1300];
% location = earthLocation(1:2);
% 
% 
% % equCoords= wrap2pi(celestial.coo.horiz_coo([az alt],julDate,location,'e'));
% % RA = equCoords(:,1); 
% % dec = equCoords(:,2);
% 
% 
% 
% [equCoords] = celestial.coo.coco([long lat],'g','j2000.0','r','r');
% RA = wrap22pi(equCoords(:,1)); %right ascension
% dec = wrap2pi(equCoords(:,2)); %declination
% 
% horzCoords= wrap2pi(celestial.coo.horiz_coo([RA dec],julDate,location,'h'));
% az = horzCoords(:,1);
% alt = horzCoords(:,2);
% 
% figure
% gridDelaunay = delaunay(az,alt);
% h = trisurf(gridDelaunay,rad2deg(az),rad2deg(alt),map.*0,map);
% set(h, 'LineStyle', 'none')
% axis equal
% % axis off
% view([0,90])
% xlabel('Az^\circ')
% ylabel('Alt^\circ')
% 
% figure
% gridDelaunayRAdec = delaunay(RA,dec);
% hRAdec = trisurf(gridDelaunayRAdec,rad2deg(RA),rad2deg(dec),map.*0,map);
% set(hRAdec, 'LineStyle', 'none')
% axis equal
% view([0,90])
% xlabel('RA^\circ')
% ylabel('Dec^\circ')
% 
% % figure
% % plot(az,alt,'.','markerSize',1)


fPlot = 408;

G0 = GlobalSkyModel2016;
G0 = G0.generate(fPlot);
G1 = G0.underSample(3);
figure,G1.plotProj(1,true,'m')
G1.plotHorizon('m')

G2 = G1.setCoorSys('RAdec');
G2 = G2.setLocation([0,0,0]);
G2 = G2.setTime(datetime(2021,09,21,0,0,0));
figure,G2.plotProj(1,true,'m')
G2.plotHorizon

G1 = G1.setTime(datetime(2021,09,21,0,0,0));
G1 = G1.setLocation([0,0,0]);
G3 = G1.setCoorSys('Horiz');
figure,G3.plotProj(1,true,'m')
G3.plotHorizon
G3.plotSun