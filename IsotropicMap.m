classdef IsotropicMap < GlobalSkyModelBase
    
    properties 

    end
    
    properties (SetAccess = private)
        spectral_index(1,1) double {mustBeReal,mustBeFinite,mustBeNegative} = -2.6
        data(:,1) double        
        Tcmb(1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 2.7
        Tg0 = 20; %K
        v0 = 408; %MHz   
    end
    
    methods 
        function obj = IsotropicMap(freq_unit,spectral_index, Tg0, v0)
            % Isotropic class constructor method

            % Inputs
            % - freq_unit: {'Hz',('MHz'),'GHz'}
            % - spectral_index: negative real number (-2.6)
            %
            % Outputs
            % - obj:    Haslam object
            %
            % Dependencies
            % -
            %
            % Created: 2021-09-06, Carla Pieterse
            % Updated: 2021-09-06, Carla Pieterse
            %
            % Tested : Matlab R2021a
            %  Level : 1
            %   File : 
            %
            % Example
            %   Hmap = Haslam;
            %   Hmap = Hmap.generate(408);
            %   Hamp.view(1,true)
        
            if nargin > 0 && ~isempty(freq_unit), obj.freq_unit = freq_unit; end
            if nargin > 1 && ~isempty(spectral_index), obj.spectral_index = spectral_index; end
            if nargin > 2 && ~isempty(spectral_index), obj.Tg0 = Tg0; end
            if nargin > 3 && ~isempty(spectral_index), obj.v0 = v0; end
            % Set the path
            %P = fileparts(mfilename('fullpath'));
            %obj.dataPath = [P,'\data\haslam408_dsds_Remazeilles2014.fits'];
            
            % Load the data
            %d_ =  fitsread(obj.dataPath,'BinaryTable');
            %d_ = transpose(d_{1});
            obj.data = ones(3145728,1);
        end
        
        function [obj, map_out] = generate(obj,freqs)
            %     [obj, map_out] = generate(obj,freqs) 
            %     Generate a global sky model at a given frequency or frequencies
            %     
            %     Parameters
            %     ----------
            %     freq: scalar or array
            %     Frequency for which to return GSM model (unit as in the object)
            %     
            %     Returns
            %     -------
            %     map_out
            %     Global sky model in healpix format, with NSIDE=512. Output map
            %     is in galactic coordinates, and in antenna temperature units (K).
            
            assert(min(size(freqs))  == 1, 'freqs must be vector')
            freqs_mhz = freqs.*obj.freqScale;
            map_out = ((obj.Tg0.*(obj.v0./freqs).^(-obj.spectral_index) + obj.Tcmb)).*obj.data;
            %map_out = (obj.data - obj.Tcmb).*(freqs_mhz./408.0).^(obj.spectral_index) + obj.Tcmb;
            
            obj.generated_map_data = map_out;
            obj.generated_map_freqs = freqs;
        end
        
        function obj = setTcmb(Tcmb)
            % SETTCMB sets the CMB temperature in the object
            % obj = setTcmb(Tcmb)
            % If Tcmb = 0 the full temperature is scled by the power law.
            % If Tcmb > 0, it is removed from the scaling and added back.
            % Default is 2.73 K
            
            obj.Tcmb = Tcmb;
            if ~isempty(obj.generated_map_data)
                obj = obj.generate(obj.generated_map_freqs);
            end
        end
    end
    
end