classdef GlobalSkyModelBase
    % Shared base superclass for GlobalSkyModel, GlobalSkyModel2016 and
    % Haslam classes
    properties
        projectionType = 'm'    % See MAAT toolbox celestial.proj.projectcoo for details (top/bot is direction cosine hemispheres)
        projectionR = 1         % See MAAT toolbox celestial.proj.projectcoo for details
        projectionPar1 = []     % See MAAT toolbox celestial.proj.projectcoo for details
        projectionPar2 = []     % See MAAT toolbox celestial.proj.projectcoo for details
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
    end
    
    properties (Dependent = true)
        Nside 
        Nf
        xyHorizon
        xySun
        julDate
    end
    
    properties (Constant = true, Hidden = true)
        signPhi = -1;
        astroGrids = {'ENU','Horiz','RAdec','GalLongLat'}
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
            Nside = sqrt(size(obj.generated_map_data,1)./12);
        end
        
        function Nf = get.Nf(obj)
            Nf = size(obj.generated_map_data,2);
        end
        
        function longlat = get.longlat(obj)
            tp = pix2ang(obj.Nside);
            tp = [tp{:}];
            th = tp(1,:);
            ph = obj.signPhi.*tp(2,:);
            longlat = [wrap2pi(ph(:)),pi/2 - th(:)];
        end
        
        function xyHorizon = get.xyHorizon(obj)
            az = linspace(-pi,pi,1001).';
            alt = zeros(size(az));
            if strcmp(obj.gridType,'Horiz')
                xyHorizon = [az,alt];
            else
                xyHorizon = wrap2pi(celestial.coo.horiz_coo([az,alt],obj.julDate,deg2rad(fliplr(obj.location(1:2))),'e'));
                if strcmp(obj.gridType,'GalLongLat')
                    xyHorizon = celestial.coo.coco(xyHorizon,'j2000.0','g','r','r');
                end
            end
            xyHorizon = wrap2pi(xyHorizon);
        end
        
        function xySun = get.xySun(obj)
            sunStruct = celestial.SolarSys.get_sun(obj.julDate,deg2rad(fliplr(obj.location(1:2))));
            xySun = [sunStruct.RA,sunStruct.Dec];
            switch obj.gridType
                case 'Horiz'
                    xySun = [sunStruct.Az,sunStruct.Alt];
                case 'GalLongLat'
                    xySun = celestial.coo.coco(xySun,'j2000.0','g','r','r');
            end
            xySun = wrap2pi(xySun);
            
        end
        
        function julDate = get.julDate(obj)
            julDate = convert.date2jd([obj.UTCtime.Day,obj.UTCtime.Month,obj.UTCtime.Year,obj.UTCtime.Hour,obj.UTCtime.Minute,obj.UTCtime.Second]);
        end
        
        function obj = setTime(obj,UTCtime)
            % SETTIME sets the time
            % obj = setTime(obj,UTCtime)
            
            obj.UTCtime = UTCtime;
            obj = obj.changeGrid(obj.gridType);  % This always goes from gal -> whereever, so update will happen automatically
        end
        
        function obj = setLocation(obj,location)
            % SETLOCATION sets the location of the observer
            % obj = = setLocation(obj,location)
            
            obj.location  = location;
            obj = obj.changeGrid(obj.gridType);  % This always goes from gal -> whereever, so update will happen automatically
        end
        
        function obj = changeGrid(obj,gridType)
            % CHANGEGRID sets the coordinate system grid
            % obj = changeGrid(obj,gridType)
            % Input: changeGrid can be anything in obj.astroGrids
            
            assert(ismember(gridType,obj.astroGrids),'Unkown coorSys. See obj.astroGrids for allowable names')
            
            if strcmp(gridType,'GalLongLat')
                obj.xy = obj.longlat;
            else %if strcmp(coorSys,'RAdec') || strcmp(coorSys,'Horiz') 
                % Always calculate this - needed for both transforms
                equCoords = celestial.coo.coco([obj.longlat],'g','j2000.0','r','r');
                obj.xy = [wrap2pi(equCoords(:,1)),wrap2pi(equCoords(:,2))];
                if strcmp(gridType,'Horiz')  % Update if needed
                    horzCoords = wrap2pi(celestial.coo.horiz_coo([obj.xy],obj.julDate,deg2rad(fliplr(obj.location(1:2))),'h'));
                    obj.xy = horzCoords(:,1:2);
                end
            end
            obj.gridType = gridType;
        end
        
        function obj = underSample(obj,sampleFactor)
           % UNDERSAMPLE returns the objects on a coarses grid
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
        end
        
        function [xy,z] = project(obj,longlat)
            % PROJECT projects the spherical map to 2D plane
            % [xy] = project(obj,longlat)
            % Uses MAAT toolbox celestial.proj.projectcoo
            % Input  : - Matrix of [Longitude,Latitude], in radians.
            % Output : - Matrix of [x,y] positions
            
            assert(size(longlat,2)==2,'Expect 2 column matrix for longlat')
            
            if any(strcmp(obj.projectionType,{'top','bot'}))
                [x,y,z] = PhTh2DirCos(-longlat(:,1),pi/2 - longlat(:,2));
                xy = [x(:),y(:)];
            else
                if contains('amhpslbtPGBCxr',obj.projectionType) || isempty(obj.projectPar1)
                    Nin = 4;
                elseif contains('cgMoS',obj.projectionType) || isempty(obj.projectPar2)
                    Nin = 5;
                else
                    Nin = 6;
                end
                switch Nin
                    case 4
                        [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),obj.projectionType,obj.projectionR);
                    case 5
                        [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),obj.projectionType,obj.projectionR,par1);
                    case 6
                        [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),obj.projectionType,obj.projectionR,obj.projectionPar1,obj.projectionPar2);
                end
                z = [];
            end
            
        end
        
        function plotProj(obj,idx,logged)
            % PLOTPROJ plots a projection of the map
            % plotProj(obj,idx,longlat)
            % idx is the frequency index to plot
            % See GlobalSkyModelBase.project for details on the input
            % structure
            
            if nargin < 2 || isempty(idx), idx = 1; end
            if nargin < 3 || isempty(logged), logged = false; end
            
            
            gmap = obj.generated_map_data(:,idx);
            if logged, gmap = log2(gmap); end
            
            if isempty(obj.xy), obj.xy = obj.longlat; end
            
            [xyProj,zProj] = obj.project(obj.xy);
            if strcmp(obj.projectionType,'top')
                xyProj = xyProj(zProj >= 0,:);
                gmap = gmap(zProj >= 0);
            elseif strcmp(obj.projectionType,'bot')
                xyProj = xyProj(zProj <= 0,:);
                gmap = gmap(zProj <= 0);
            end
%             plot3(xyProj(:,1),xyProj(:,2),gmap,'.')
            gridDelaunay = delaunay(xyProj(:,1),xyProj(:,2));
            h = trisurf(gridDelaunay,xyProj(:,1),xyProj(:,2),gmap.*0,gmap);
            set(h, 'LineStyle', 'none')
            axis equal
            % axis off
            view([0,90])
            hold on
        end
        
        function plotHorizon(obj,style)
            % PLOTHORIZON plots the horizon on the current axis
            % plotHorizon(obj,style)
            % style is the linestyle of the plot ('w-')
            
            if nargin < 2 || isempty(style), style = 'w-'; end
            
            xyProj = obj.project(obj.xyHorizon);
            plot(xyProj(:,1),xyProj(:,2),style)
        end
        
        function plotSun(obj,style)
            % PLOTSUN plots the approximate sun position
            % plotSun(obj,style)
            % style is the markerstyle of the plot ('w*')
            
            if nargin < 2 || isempty(style), style = 'w*'; end
            
            [xyProj,zProj] = obj.project(obj.xySun);
            if strcmp(obj.projectionType,'top')
                xyProj = xyProj(zProj >= 0,:);
            elseif strcmp(obj.projectionType,'bot')
                xyProj = xyProj(zProj <= 0,:);
            end
            plot(xyProj(:,1),xyProj(:,2),style)

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