classdef GlobalSkyModelBase
    % Shared base superclass for GlobalSkyModel and GlobalSkyModel2016
    properties (SetAccess = protected)
        freq_unit(1,:) char {mustBeMember(freq_unit,{'Hz','kHz','MHz','GHz'})} = 'MHz'
        generated_map_data(:,:) double = []
        generated_map_freqs(1,:) double = []
    end
    
    properties (SetAccess = protected, Hidden = true)
        dataPath
    end
         
    properties (Dependent = true, Hidden = true)
       freqScale 
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