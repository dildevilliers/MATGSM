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
    
    properties (Constant = true, Hidden = true)
        kB = 1.38065e-23
        C = 2.99792e8
        h = 6.62607e-34
        T = 2.725
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
            
            % Take transposes below - to be the same as the pygsm python which
            % reads in a transposed matrix from h5 files.
            map_ni_ = obj.map_ni.';
            spec_nf_ = obj.spec_nf.';
            
            nfreq = size(spec_nf_,2);
            
            map_out = zeros(numel(freqs_ghz), size(map_ni_,2));
            for ff = 1:numel(freqs_ghz)
                freq = freqs_ghz(ff);
                left_index = -1;
                for ii = 1:(nfreq-1)
                    if spec_nf_(1,ii) <= freq && freq <= spec_nf_(1,ii+1)
                        left_index = ii;
                        break;
                    end
                end
                
                interp_spec_nf = spec_nf_;
                interp_spec_nf(1:2,:) = real(log10(interp_spec_nf(1:2,:)));  % This only goes to index 2 - the python in pygsm2016 want to go to 3 but that can't be (negative values in the argument for the third row...)
                
                % Don't worry about the left_index = -1 possibility. It should never get to
                % that due to the earlier asserts forcing the evaluated
                % frequency into the right ranges
                x1 = interp_spec_nf(1,left_index);
                x2 = interp_spec_nf(1,left_index+1);
                y1 = interp_spec_nf(2:end,left_index);
                y2 = interp_spec_nf(2:end,left_index+1);
                x = log10(freq);
                interpolated_vals = (x.*(y2 - y1) + x2.*y1 - x1.*y2)./(x2 - x1);
                map_out(ff,:) = sum(bsxfun(@times,10.^(interpolated_vals(1)).*interpolated_vals(2:end),map_ni_),1);
                
                nMap = size(map_out,2);
                nside = sqrt(nMap./12);
                if ff == 1, iRing = ring2nest(nside,1:nMap); end
                map_out(ff,:) = map_out(ff,iRing);
                
                switch obj.unit
                    case 'TCMB'
                        conversion = 1./obj.K_CMB2MJysr(1,1e9.*freq);
                    case 'TRJ'
                        conversion = 1./obj.K_RJ2MJysr(1,1e9.*freq);
                    otherwise
                        conversion = 1;
                end
                map_out(ff,:) = map_out(ff,:).*conversion;
            end
            
            obj.generated_map_freqs = freqs;
            obj.generated_map_data = map_out.';
        end
        
        % Conversion factor functions
        
        function MJysr = K_CMB2MJysr(obj,K_CMB, nu) % in Kelvin and Hz
            hoverk = obj.h./obj.kB;
            B_nu = 2.*(obj.h.*nu).*(nu./obj.C).^2./(exp(hoverk.*nu./obj.T) - 1);
            conversion_factor = (B_nu.*obj.C./nu./obj.T).^2./2.*exp(hoverk.*nu./obj.T)./obj.kB;
            MJysr = K_CMB.*conversion_factor.*1e20; %1e-26 for Jy and 1e6 for MJy
        end
        
        function MJysr = K_RJ2MJysr(obj,K_RJ, nu) %in Kelvin and Hz
            conversion_factor = 2.*(nu./obj.C).^2.*obj.kB;
            MJysr = K_RJ.*conversion_factor.*1e20; %1e-26 for Jy and 1e6 for MJy
        end
    end
end


