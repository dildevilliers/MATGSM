close all
clear all

GSM = GlobalSkyModel([],'haslam');
GSM = GSM.generate([408]);
GSM.view(1,true)
GSM.write_fits('c:\Temp\GSMfits.fits')

%%
signPhi = -1;
map = GSM.generated_map_data(:,1);
sz = size(map);
nside = sqrt(max(sz)/12);
tp = pix2ang(nside);
tp = [tp{:}];
th = tp(1,:);
ph = signPhi.*tp(2,:);
long = pi/2 - th(:);
lat = wrap2pi(ph(:));

time = datetime(2018,7,22,0,0,0);
julDate = convert.date2jd([time.Day,time.Month,time.Year,time.Hour,time.Minute,time.Second]);
earthLocation = [deg2rad(18.86) deg2rad(-33.93) 300];
location = earthLocation(1:2);


% equCoords= wrap2pi(celestial.coo.horiz_coo([az alt],julDate,location,'e'));
% RA = equCoords(:,1); 
% dec = equCoords(:,2);



[equCoords] = celestial.coo.coco([long lat],'g','j2000.0','r','r');
RA = wrap2pi(equCoords(:,1)); %right ascension
dec = wrap2pi(equCoords(:,2)); %declination

horzCoords= wrap2pi(celestial.coo.horiz_coo([RA dec],julDate,location,'h'));
az = horzCoords(:,1);
alt = horzCoords(:,2);

figure
gridDelaunay = delaunay(az,alt);
h = trisurf(gridDelaunay,az,alt,map.*0,map);

set(h, 'LineStyle', 'none')
axis equal
axis off
view([0,90])

figure
plot(az,alt,'.','markerSize',1)

