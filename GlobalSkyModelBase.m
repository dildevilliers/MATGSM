classdef GlobalSkyModelBase
    % Shared base superclass for GlobalSkyModel, GlobalSkyModel2016 and
    % Haslam classes
    properties (SetAccess = protected)
        freq_unit(1,:) char {mustBeMember(freq_unit,{'Hz','kHz','MHz','GHz'})} = 'MHz'
        generated_map_data(:,:) double = []     % Always in galactic coordinates, on a HEALpix ring sampling grid
        generated_map_freqs(1,:) double = []
    end
    
    properties (SetAccess = private)
        coorSys = 'galactic'        % Can be set to {galactic,equatorial,horizontal} with projector
    end
    
    properties (SetAccess = protected, Hidden = true)
        dataPath
    end
         
    properties (Dependent = true, Hidden = true)
       freqScale 
       longlat  % Longitude and latitude of the grid (N x 2)
    end
    
    properties (Dependent = true)
        Nside 
        Nf
    end
    
    properties (Constant = true, Hidden = true)
        signPhi = -1;
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
        
        function [RAdec] = coor2RAdec(obj)
            % COOR2RADEC returns the [RA,Dec] coordinates of the object
            equCoords = celestial.coo.coco([obj.longlat],'g','j2000.0','r','r');
            RAdec = [wrap22pi(equCoords(:,1)),wrap2pi(equCoords(:,2))]; 
        end
        
        function [azAlt,UTCtime,location] = coor2hor(obj,UTCtime,location)
            % COOR2HOR returns the horizontal (az,alt) coordinates of the
            % object
            % [azAlt,locTime,location] = coor2hor(obj,UTCtime,location)
            %
            % Inputs:
            % UTCtime in datetime format (default = 2018/01/01 , 00:00)
            % location as [long (rad) lat (rad) alt (m)] with default = [0.3292   -0.5922  300.0000] 
            
            if nargin < 2 || isempty(UTCtime), UTCtime = datetime(2018,1,1,0,0,0); end
            if nargin < 3 || isempty(location), location = [0.3292   -0.5922  300.0000]; end
            
            julDate = convert.date2jd([UTCtime.Day,UTCtime.Month,UTCtime.Year,UTCtime.Hour,UTCtime.Minute,UTCtime.Second]);
            horzCoords = wrap2pi(celestial.coo.horiz_coo([obj.coor2RAdec],julDate,location(1:2),'h'));
            azAlt = horzCoords(:,1:2);
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
end