classdef IsotropicMap < GlobalSkyModelBase
    
    properties 

    end
    
    properties (SetAccess = private)
        spectral_index(1,1) double {mustBeReal,mustBeFinite,mustBeNegative} = -2.6
        data(:,1) double        
        Tcmb(1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 2.73
        Tg0(1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 20; % K
        v0(1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 408; % MHz   
    end

    properties (Dependent = true)
        Npix
    end
    
    methods 
        function obj = IsotropicMap(freq_unit, spectral_index, Tg0, v0, resIdx)
            % Isotropic class constructor method

            % Inputs
            % - freq_unit: {'Hz',('MHz'),'GHz'}
            % - spectral_index: negative real number (-2.6)
            % - Tg0: Reference temperature in (K) at...
            % - v0: the reference frequency in MHz
            % - resIdx: the resolution index (integer between 1 and 10)
            %
            % Outputs
            % - obj:    IsotropicMap object
            %
            % Dependencies
            % -
            %
            % Created: 2021-09-06, Carla Pieterse
            % Updated: 2023-06-08, Dirk de  Villiers
            %
            % Tested : Matlab R2021a
            %  Level : 1
            %   File : 
            %
            % Example
            %   Imap = IsotropicMap;
            %   Imap = Imap.generate(408);
            %   Iamp.view(1,true)
        
            if nargin > 0 && ~isempty(freq_unit), obj.freq_unit = freq_unit; end
            if nargin > 1 && ~isempty(spectral_index), obj.spectral_index = spectral_index; end
            if nargin > 2 && ~isempty(Tg0), obj.Tg0 = Tg0; end
            if nargin > 3 && ~isempty(v0), obj.v0 = v0; end
            if nargin > 4 && ~isempty(resIdx)
                assert(mod(resIdx,1) == 0 && resIdx > 0,'resIdx must be a positive integer')
                rI = resIdx;
            else
                rI = 8;
            end
            Npixels = 12.*(2.^rI).^2;
            
            obj.data = ones(Npixels,1);
        end

        function Npix = get.Npix(obj)
            if isempty(obj.generated_map_data)
                Npix = size(obj.data,1);
            else
                Npix = size(obj.generated_map_data,1);
            end
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
            freqs = freqs(:).';   % Force row vector
            freqs_mhz = freqs.*obj.freqScale;
            map_out = ((obj.Tg0.*(freqs_mhz./obj.v0).^(obj.spectral_index) + obj.Tcmb)).*obj.data;
            
            obj.generated_map_data = map_out;
            obj.generated_map_freqs = freqs;
        end
        
        function obj = setTcmb(obj,Tcmb)
            % SETTCMB sets the CMB temperature in the object
            % obj = setTcmb(obj,Tcmb)
            % If Tcmb = 0 the full temperature is scaled by the power law.
            % If Tcmb > 0, it is removed from the scaling and added back.
            % Default is 2.73 K
            
            obj.Tcmb = Tcmb;
            if ~isempty(obj.generated_map_data)
                obj = obj.generate(obj.generated_map_freqs);
            end
        end

    end
    
end