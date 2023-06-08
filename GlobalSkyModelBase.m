classdef GlobalSkyModelBase
    % Shared base superclass for GlobalSkyModel, GlobalSkyModel2016 and
    % Haslam classes
    properties
        projectionType = 'm'    % See MAAT toolbox celestial.proj.projectcoo for details 
        projectionR = 1         % See MAAT toolbox celestial.proj.projectcoo for details
        projectionPar1 = []     % See MAAT toolbox celestial.proj.projectcoo for details
        projectionPar2 = []     % See MAAT toolbox celestial.proj.projectcoo for details
        
        hideGround(1,1) logical = false     % Set to true to hide the ground region in all plots
        gridPlotStepDeg(1,1) double = 10;   % Density of grid lines for plots
    end
    
    properties (SetAccess = protected)
        freq_unit(1,:) char {mustBeMember(freq_unit,{'Hz','kHz','MHz','GHz'})} = 'MHz'
        generated_map_data(:,:) double = []     % Always in galactic coordinates, on a HEALpix ring sampling grid
        generated_map_freqs(1,:) double = []
    end
    
    properties (SetAccess = private)
        gridType = 'GalLongLat'        % Can be set to {'Horiz','RAdec','GalLongLat'} setCoorSys
        location(1,3) double = [(-30.721745), (21.411701),  300.0000]  % Earth location in [Lat(deg) Long(deg) mASL]
        UTCtime(1,1) datetime = datetime(2018,1,1,0,0,0)
    end
    
    properties (SetAccess = protected, Hidden = true)
        dataPath
    end
         
    properties (SetAccess = private, Hidden = true)
        xy(:,2) double  % The current local grid
    end
    
    properties (Dependent = true, Hidden = true)
       freqScale 
       longlat  % Longitude and latitude of the original galactic grid (N x 2)
       xyHorizon
       xySunMoon
       xyGridAz 
       xyVerifyMarkers
       idxGround
    end
    
    properties (Dependent = true)
        Nside 
        resIdx
        Nf
        signPhi
        julDate
    end
    
    properties (Constant = true, Hidden = true)
        astroGrids = {'Horiz','RAdec','GalLongLat'}
        verifyMarkers = {'GC','VelaSNR','Cygnus','Cas-A','Cen-A','Tau-A','Orion-A','LMC','SMC'}
        verifyMarkerGalCoors = [0,0; -96,-3.3; 77,2; 111.75,-2.11; -50.5,19.42; -175.42,-5.79; -151,-19.36; -79.5,-32.85; -57.2,-44.3];
    end
    
    methods
        function freqScale = get.freqScale(obj)
            switch obj.freq_unit
                case 'Hz'
                    freqScale = 1e-6;
                case 'kHz'
                    freqScale = 1e-3;
                case 'MHz'
                    freqScale = 1;
                case 'GHz'
                    freqScale = 1e3;
            end
        end
        
        function Nside = get.Nside(obj)
            Nside = nPix2nSide(obj.Npix);
        end

        function resIdx = get.resIdx(obj)
            resIdx = log2(obj.Nside);
        end

        function Nf = get.Nf(obj)
            Nf = size(obj.generated_map_data,2);
        end
        
        function signPhi = get.signPhi(obj)
            switch obj.gridType
                case 'GalLongLat'
                    signPhi = -1;
                case 'RAdec'
                    signPhi = -1;
                case 'Horiz'
                    signPhi = 1;
            end
        end
        
        function longlat = get.longlat(obj)
            tp = pix2ang(obj.Nside);
            tp = [tp{:}];
            th = tp(1,:);
            ph = tp(2,:);
            longlat = [wrap2pi(ph(:)),pi/2 - th(:)];
        end
        
        function idxGround = get.idxGround(obj)
            obj = obj.changeGrid('Horiz');
            idxGround = find(obj.xy(:,2) < 0);
        end
        
        function xyHorizon = get.xyHorizon(obj)
            % Work in horizontal coordinates
            % First four elements are the cardinal points for
            % plotDirections
            az = [deg2rad(0:90:270).';linspace(-pi,pi,1001).'];
            alt = zeros(size(az));
            
            xyHorizon = [az,alt];
            if ~strcmp(obj.gridType,'Horiz')
                % Go to equatorial
                xyHorizon = wrap2pi(horiz_coo(xyHorizon,obj.julDate,deg2rad(fliplr(obj.location(1:2))),'e'));
                if strcmp(obj.gridType,'GalLongLat')
                    % Go to galactic
                    xyHorizon = celestial.coo.coco(xyHorizon,'j2000.0','g','r','r');
                end
            end
            xyHorizon = wrap2pi([obj.signPhi,1].*xyHorizon);
        end
        
        function xyGridAz = get.xyGridAz(obj)
            % Make vectors of required grid points
            xyGridAz.latVect = deg2rad(-90+obj.gridPlotStepDeg:obj.gridPlotStepDeg:90-obj.gridPlotStepDeg);
            xyGridAz.longVect = deg2rad(-180+obj.gridPlotStepDeg:obj.gridPlotStepDeg:180);
            xyGridAz.Nlat = length(xyGridAz.latVect);
            xyGridAz.Nlong = length(xyGridAz.longVect);
            xyGridAz.Nlines = xyGridAz.Nlat + xyGridAz.Nlong;
            
            % Pack the lines as increasing parallels, and then increasing
            % meridians
            Nplot = 201;
            gridLines = zeros(Nplot,2,xyGridAz.Nlines);
            for nnLat = 1:xyGridAz.Nlat
                az = linspace(-pi,pi,Nplot).';
                alt = ones(size(az)).*xyGridAz.latVect(nnLat);
                gridLines(:,:,nnLat) = [az,alt];
            end
            for nnLong = 1:xyGridAz.Nlong
                alt = linspace(-pi/2,pi/2,Nplot).';
                az = ones(size(alt)).*xyGridAz.longVect(nnLong);
                gridLines(:,:,xyGridAz.Nlat + nnLong) = [az,alt];
            end
           
            gridLines = reshape(gridLines,Nplot*xyGridAz.Nlines,2);
            if ~strcmp(obj.gridType,'Horiz')
                gridLines = wrap2pi(horiz_coo(gridLines,obj.julDate,deg2rad(fliplr(obj.location(1:2))),'e'));
                if strcmp(obj.gridType,'GalLongLat')
                    gridLines = celestial.coo.coco(gridLines,'j2000.0','g','r','r');
                end
            end
            gridLines = reshape(gridLines,Nplot,2,xyGridAz.Nlines);
            xyGridAz.gridLines = wrap2pi([obj.signPhi,1].*gridLines);
        end
        
        function xySunMoon = get.xySunMoon(obj)
            sunStruct = celestial.SolarSys.get_sun(obj.julDate,deg2rad(fliplr(obj.location(1:2))));
            moonStruct = celestial.SolarSys.get_moon(obj.julDate,deg2rad(fliplr(obj.location(1:2))));
            xySunMoon = [sunStruct.RA,sunStruct.Dec;moonStruct.RA,moonStruct.Dec];
            switch obj.gridType
                case 'Horiz'
                    xySunMoon = [sunStruct.Az,sunStruct.Alt;moonStruct.Az,moonStruct.Alt];
                case 'GalLongLat'
                    xySunMoon = celestial.coo.coco(xySunMoon,'j2000.0','g','r','r');
            end
            xySunMoon = wrap2pi([obj.signPhi,1].*xySunMoon);
        end
        
        function xyVerifyMarkers = get.xyVerifyMarkers(obj)
            
            xyVerifyMarkersGC = deg2rad(obj.verifyMarkerGalCoors);
            if strcmp(obj.gridType,'GalLongLat')
                xyVerifyMarkers = wrap2pi([obj.signPhi,1].*xyVerifyMarkersGC);
            else 
                xyVerifyMarkersEq = celestial.coo.coco(xyVerifyMarkersGC,'g','j2000.0','r','r');
                xyVerifyMarkers = wrap2pi([obj.signPhi,1].*xyVerifyMarkersEq);
                if any(strcmp(obj.gridType,{'Horiz'}))  % Update if needed
                    xyVerifyMarkersHor = wrap2pi(celestial.coo.horiz_coo(xyVerifyMarkersEq,obj.julDate,deg2rad(fliplr(obj.location(1:2))),'h'));
                    xyVerifyMarkers =  wrap2pi([obj.signPhi,1].*xyVerifyMarkersHor);
                end
            end
        end
        
        function julDate = get.julDate(obj)
            julDate = convert.date2jd([obj.UTCtime.Day,obj.UTCtime.Month,obj.UTCtime.Year,obj.UTCtime.Hour,obj.UTCtime.Minute,obj.UTCtime.Second]);
        end
        
        function obj = setTime(obj,UTCtime)
            % SETTIME sets the time
            % obj = setTime(obj,UTCtime)
            
            obj.UTCtime = UTCtime;
            if ~isempty(obj.xy)
                obj = obj.changeGrid(obj.gridType);  % This always goes from gal -> whereever, so update will happen automatically
            end
        end
        
        function obj = setLocation(obj,location)
            % SETLOCATION sets the location of the observer
            % obj = setLocation(obj,location)
            % location is a 3 element vector: [Lat(deg) Long(deg) mASL]
            
            %obj.location  = location;
            obj.location = setLocation(location);
            if ~isempty(obj.xy)
                obj = obj.changeGrid(obj.gridType);  % This always goes from gal -> whereever, so update will happen automatically
            end
        end
        
        function obj = changeGrid(obj,gridType)
            % CHANGEGRID sets the coordinate system grid
            % obj = changeGrid(obj,gridType)
            % Input: changeGrid can be anything in obj.astroGrids
            
            assert(ismember(gridType,obj.astroGrids),'Unkown coorSys. See obj.astroGrids for allowable names')
            
            obj.gridType = gridType;
            
            if strcmp(gridType,'GalLongLat')
                obj.xy = [obj.signPhi,1].*obj.longlat;
            else %if any(strcmp(coorSys,{'RAdec','Horiz'}))
                % Always calculate this - needed for all three transforms
                xyEq = celestial.coo.coco(obj.longlat,'g','j2000.0','r','r');
                obj.xy = wrap2pi([obj.signPhi,1].*xyEq);
                if any(strcmp(gridType,{'Horiz'}))  % Update if needed
                    xyHor = wrap2pi(horiz_coo(xyEq,obj.julDate,deg2rad(fliplr(obj.location(1:2))),'h'));
                    obj.xy = wrap2pi([obj.signPhi,1].*xyHor);
                end
            end
        end
        
        function obj = underSample(obj,sampleFactor)
           % UNDERSAMPLE returns the objects on a coarser grid
           % obj = underSample(obj,sampleFactor)
           % SampleFactor is a positive integer indicating 2^sampleFactor reduction
           % in nSide with default (1)
           
           if nargin < 2 || isempty(sampleFactor), sampleFactor  = 1; end
           assert(sampleFactor >= 0 && ~mod(sampleFactor,1),'Expected integer value for sampleFactor') 
           
           % Iteratively make the thing smaller
           for ii = 1:sampleFactor
               % Get nested indexes
               L = size(obj.generated_map_data,1);
               iNest = nest2ring(obj.Nside,(1:L).');
               dNest = zeros(L/4,obj.Nf);
               for ff = 1:obj.Nf
                   % Taking the raw average here - reshape to 4 rows, add them
                   % and normalise by 4
                   dNest(:,ff) = (sum(reshape(obj.generated_map_data(iNest,ff),4,L/4),1)./4).';
               end
               iRing = ring2nest(obj.Nside./2,(1:nSide2nPix(obj.Nside./2)).');
               obj.generated_map_data = dNest(iRing,:);
           end
           obj = obj.changeGrid(obj.gridType);  % Set xy
        end
        
        function obj = updateMapData(obj,mapData)
            % UPDATEMAPDATA updates the results in generated_map_data
            % Usually from external calculations after internally setting
            % up the result
            
            assert(all(size(mapData) == size(obj.generated_map_data)),'Input mapData must be the same size as internal generated_map_data')
            
            obj.generated_map_data = mapData;
        end
        
        function [xy,z] = project(obj,longlat)
            % PROJECT projects the spherical map to 2D plane
            % [xy] = project(obj,longlat)
            % Uses MAAT toolbox celestial.proj.projectcoo
            % Input  : - Matrix of [Longitude,Latitude], in radians.
            % Output : - Matrix of [x,y] positions
            
            assert(size(longlat,2)==2,'Expect 2 column matrix for longlat')
            
            if contains('amhpslbtPGBCxr',obj.projectionType) || isempty(obj.projectionPar1)
                Nin = 4;
            elseif contains('cgMoS',obj.projectionType) || isempty(obj.projectionPar2)
                Nin = 5;
            else
                Nin = 6;
            end
            switch Nin
                case 4
                    [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),obj.projectionType,obj.projectionR);
                case 5
                    [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),obj.projectionType,obj.projectionR,obj.projectionPar1);
                case 6
                    [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),obj.projectionType,obj.projectionR,obj.projectionPar1,obj.projectionPar2);
            end
            z = [];  % for in case later...
        end
        
        function plotProj(obj,idx,logged)
            % PLOTPROJ plots a projection of the map
            % plotProj(obj,idx,logged)
            % idx is the frequency index to plot
            % logged is logical to plot in log scale (true)
            
            if nargin < 2 || isempty(idx), idx = 1; end
            if nargin < 3 || isempty(logged), logged = true; end
            
            
            gmap = obj.generated_map_data(:,idx);
            if logged, gmap = log2(gmap); end
            
            if obj.hideGround, gmap(obj.idxGround) = nan; end
            
            if isempty(obj.xy)  % we are in Galactic coordinates now
                obj.xy = [obj.signPhi,1].*obj.longlat; 
            end
            
            [xyProj] = obj.project(obj.xy);
            gridDelaunay = delaunay(xyProj(:,1),xyProj(:,2));
            h = trisurf(gridDelaunay,xyProj(:,1),xyProj(:,2),gmap.*0,gmap);
            set(h, 'LineStyle', 'none')
            axis equal
            axis off
            view([0,90])
            hold on
            colorbar('location','SouthOutside')
        end
        
        function plotSkyView(obj,idx,logged,details)
            % PLOTSKYVIEW plots the sky view at the current time and place
            % plotSkyView(obj,idx,logged,paramsPlot)
            % idx is the frequency index to plot
            % logged is logical to plot in log scale (true)
            % details is logicals requesting plotting of 
            %   [Directions, Sun, Moon, verifyMarkers]
            
            if nargin < 2 || isempty(idx), idx = 1; end
            if nargin < 3 || isempty(logged), logged = true; end
            if nargin < 4 || isempty(details), details = [0,0,0,0]; end
            
            if numel(details) < 4
                details_ = zeros(1,4);
                details_(1:numel(details)) = details;
                details = details_;
            end
            
            gmap = obj.generated_map_data(:,idx);
            if logged, gmap = log2(gmap); end
            
            if isempty(obj.xy), obj.xy = obj.longlat; end
            
            obj.projectionType = 'S';
            obj.projectionPar1 = deg2rad([0,90]);
            obj = obj.changeGrid('Horiz');
            obj.hideGround = true;
            valInd = find(obj.xy(:,2) >= 0);
            [xyProj,~] = obj.project(obj.xy(valInd,:));
            
            gridDelaunay = delaunay(xyProj(:,1),xyProj(:,2));
            h = trisurf(gridDelaunay,xyProj(:,1),xyProj(:,2),gmap(valInd).*0,gmap(valInd));
            set(h, 'LineStyle', 'none')
            axis equal
            axis off
            view([0,90])
            hold on
            colorbar('location','SouthOutside')
            if details(1), obj.plotDirections('k'); end
            if details(2) && obj.xySunMoon(1,2) >= 0, obj.plotSun; end
            if details(3) && obj.xySunMoon(2,2) >= 0, obj.plotMoon; end
            if details(4), obj.plotVerifyMarkers; end
        end
        
        function plotHorizon(obj,style)
            % PLOTHORIZON plots the horizon on the current axis
            % plotHorizon(obj,style)
            % style is the linestyle of the plot ('w-')
            
            if nargin < 2 || isempty(style), style = 'w.'; end
            
            xyProj = obj.project(obj.xyHorizon);
            plot(xyProj(5:end,1),xyProj(5:end,2),style)
        end
        
        function plotGridAz(obj,style)
            % PLOTGRIDAZ plots the azimuthal grid on the current axis
            % plotGridAz(obj,style)
            % style is the linestyle of the plot ('w-')
            % Not very reliable yet
            
            if nargin < 2 || isempty(style), style = 'w-'; end
            
            for ll = 1:obj.xyGridAz.Nlines
                xyProj = obj.project(obj.xyGridAz.gridLines(:,:,ll));
                plot(xyProj(:,1),xyProj(:,2),style,'linewidth',0.5), hold on
            end
        end
        
        function plotDirections(obj,color)
            % PLOTDIRECTIONS prints the cardinal directions on the current axis
            % plotDirections(obj,color)
            % color is the color of the text ('w')
            
            if nargin < 2 || isempty(color), color = 'w'; end
            
            textStr = 'NESW';
            xyProj = obj.project(obj.xyHorizon);
            for tt = 1:4
                text(xyProj(tt,1),xyProj(tt,2),textStr(tt),'Color',color);
            end
        end
        
        function plotSun(obj,style)
            % PLOTSUN plots the approximate sun position
            % plotSun(obj,style)
            % style is the markerstyle of the plot ('w*')
            
            if nargin < 2 || isempty(style), style = 'w*'; end
            
            sunStruct = celestial.SolarSys.get_sun(obj.julDate,deg2rad(fliplr(obj.location(1:2))));
            if ~(obj.hideGround && sunStruct.Alt < 0)  
                [xyProj,~] = obj.project(obj.xySunMoon);
                plot(xyProj(1,1),xyProj(1,2),style)
            end
        end

        function plotMoon(obj,style)
            % PLOTMOON plots the approximate moon position
            % plotMoon(obj,style)
            % style is the markerstyle of the plot ('w*')
            
            if nargin < 2 || isempty(style), style = 'wo'; end
            
            moonStruct = celestial.SolarSys.get_moon(obj.julDate,deg2rad(fliplr(obj.location(1:2))));
            if ~(obj.hideGround && moonStruct.Alt < 0)
                [xyProj,~] = obj.project(obj.xySunMoon);
                plot(xyProj(2,1),xyProj(2,2),style)
            end
        end
        
        function plotVerifyMarkers(obj,color)
            % PLOTVERIFYMARKERS plots the verification markers
            % plotVerifyMarkers(obj,color,idxPlot)
            % color is the text color ('w')
            
            if nargin < 2 || isempty(color), color = 'w'; end
            
            xyProj = obj.project(obj.xyVerifyMarkers);
            
            if obj.hideGround
                obj.gridType = 'Horiz';
                idxPlot = find(obj.xyVerifyMarkers(:,2) >= 0);
            else
                idxPlot = 1:numel(obj.verifyMarkers);
            end
            
            for tt = 1:length(idxPlot)
                text(xyProj(idxPlot(tt),1),xyProj(idxPlot(tt),2),obj.verifyMarkers(idxPlot(tt)),'Color',color);
            end
        end
            
        function view(obj, idx, logged)
            %     View generated map using mollweide projection.
            % 
            %     Parameters
            %     ----------
            %     idx: int (1)
            %         index of map to view. Only required if you generated maps at
            %         multiple frequencies.
            %     logged: logical (false)
            %         Take the log of the data before plotting. Defaults to
            %         False..
            
            
           assert(~isempty(obj.generated_map_data),'No GSM map has been generated yet. Run generate() first.')
           if nargin > 1 && ~isempty(idx)
               assert(numel(idx) == 1,'Scalar idx expected')
               gmap = obj.generated_map_data(:,idx);
               freq = obj.generated_map_freqs(idx);
           else
               gmap = obj.generated_map_data(:,1);
               freq = obj.generated_map_freqs(1);
           end
           
           if nargin < 3, logged = false; end
           if logged, gmap = log2(gmap); end
           
           healpixPlotMollweide(gmap)
%            title(['Global Sky Model at ', num2str(freq), ' MHz from the ', obj.basemap, ' map'])
        end
        
        function write_fits(obj,filename)
           fitswrite(obj.generated_map_data,filename); 
        end
    end
    
    methods (Abstract = true)
        [obj,map_out] = generate(obj,freqs)
        
    end
    
    methods (Static = true)
%         
    end
end