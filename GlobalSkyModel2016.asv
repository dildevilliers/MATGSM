classdef GlobalSkyModel2016 < GlobalSkyModelBase
    % Global sky model (GSM) class for generating sky models.
    
    properties (SetAccess = private)
        unit(1,:) char {mustBeMember(unit,{'MJysr','TCMB','TRJ'})} = 'TCMB'
        resolution(1,:) char {mustBeMember(resolution,{'hi','low'})} = 'hi'
        theta_rot(1,1) double {mustBeReal,mustBeFinite} = 0
        phi_rot(1,1) double {mustBeReal,mustBeFinite} = 0
        
        map_ni(:,:)
        spec_nf(:,:)
    end
    
    properties (SetAccess = private, Hidden = true)
        
    end
    
    methods
        function obj = GlobalSkyModel2016(freq_unit,unit,resolution,theta_rot,phi_rot)
            % GLOBALSKYMODEL2016 class constructor method
            % obj = GlobalSkyModel2016(freq_unit,unit,resolution,theta_rot,phi_rot)
            %
            % Notes
            %     -----
            %     Upon initialization, the map PCA data are loaded into memory and interpolation
            %     functions are pre-computed.
            % 
            % Inputs
            % - freq_unit: {'Hz',('MHz'),'GHz'}
            % - unit:      {'MJysr', ('TCMB'), 'TRJ'}
            % - resolution:{('hi'),'low'}
            %             Resolution of output map. Either 300 arcmin (low) or 24 arcmin (hi).
            %             For frequencies under 10 GHz, output is 48 arcmin.
            % - theta_rot: theta rotation angle in rad
            % - phi_rot: phi rotation angle in rad
            %
            % Outputs
            % - obj:    GlobalSkyModel2016 object
            %
            % Dependencies
            % -
            %
            % Created: 2021-06-10, Dirk de Villiers
            % Updated: 2021-06-10, Dirk de Villiers
            %
            % Tested : Matlab R2020a
            %  Level : 2
            %   File : testScript_GSM.m
            %
            % Example
            %   GSM = GlobalSkyModel2016;
            %   GSM = GSM.generate(408);
            %   GSM.view(1,true)
            
            if nargin > 0 && ~isempty(freq_unit), obj.freq_unit = freq_unit; end
            if nargin > 1 && ~isempty(unit), obj.unit = unit; end
            if nargin > 2 && ~isempty(resolution), obj.resolution = lower(resolution); end
            if nargin > 3 && ~isempty(theta_rot), obj.theta_rot = theta_rot; end
            if nargin > 4 && ~isempty(phi_rot), obj.phi_rot = phi_rot; end
            
            if obj.theta_rot ~= 0 || obj.phi_rot ~= 0, error('Rotation not implimented yet - set theta_rot and phi_rot to zero'); end
            
            % Set the path
            P = fileparts(mfilename('fullpath'));
            obj.dataPath = [P,'\data\gsm2016_components.h5'];
            
            % Map data to load
            labels = {'Synchrotron', 'CMB', 'HI', 'Dust1', 'Dust2', 'Free-Free'};
            fileInfo = h5info(obj.dataPath);
            switch obj.resolution
                case 'hi'
                    idxHighRes = find(strncmp({fileInfo.Datasets.Name},'highres',6));
                    obj.map_ni = zeros(fileInfo.Datasets(idxHighRes(1)).Dataspace.Size,numel(labels));
                    for ll = 1:numel(labels)
                        obj.map_ni(:,ll) = h5read(obj.dataPath,['/highres_',labels{ll},'_map']);
                    end
                case 'low'
                    obj.map_ni = h5read(obj.dataPath,'/lowres_maps');
            end
            obj.spec_nf = h5read(obj.dataPath,'/spectra');
            
            % TODO: Rotate the map...
            
            
        end
        
        function [obj, map_out] = generate(obj,freqs)
        %     Generate a global sky model at a given frequency or frequencies
        % 
        %     Parameters
        %     ----------
        %     freqs: scalar or array
        %         Frequency for which to return GSM model
        % 
        %     Returns
        %     -------
        %     map_out
        %         Global sky model in healpix format, with NSIDE=1024. Output map
        %         is in galactic coordinates, ring format.
        
        assert(min(size(freqs))  == 1, 'freqs must be vector')
        freqs_ghz = freqs.*obj.freqScale./1e3;
        assert(min(freqs_ghz) >= 0.01 && min(freqs_ghz) <= 5000, 'Frequency values lie outside 10 MHz < f < 5 THz')
        
        nfreq = size(obj.spec_nf,2);
        
        map_out = zeros(numel(freqs_ghz), size(obj.map_ni,2));
        
        for ff = 1:numel(freq_ghz)
            left_index = -1;
            for ii = 1:nfreq
                if obj.spec_nf(1,ii) <= freqs_ghz(ff) && freqs_ghz(ff) <= pbj.spec_nf(1,ii+1)
                    left_index = ii;
                    break;
                end
            end
        end
        
        keyboard

        end
    end
end


