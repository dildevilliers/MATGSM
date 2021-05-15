classdef GlobalSkyModel
    % Global sky model (GSM) class for generating sky models
    
    properties (SetAccess = private)
        freq_unit(1,:) char {mustBeMember(freq_unit,{'Hz','kHz','MHz','GHz'})} = 'MHz'
        basemap(1,:) char {mustBeMember(basemap,{'haslam','wmap','5deg'})} = '5deg'
        interpolation(1,:) char {mustBeMember(interpolation,{'cubic','pchip'})} = 'pchip'
        
        pca_map_data(:,:)
        interp_comps(1,4) cell
        generated_map_data(:,:) double = []
        generated_map_freqs(1,:) double = []
    end
    
    properties (SetAccess = private, Hidden = true)
        dataPath
        pca_freqs_mhz
        pca_scaling
        pca_comps
        generated_interp_comps(1,4) cell
    end
    
    properties (Dependent = true, Hidden = true)
       freqScale 
    end
    
    methods
        function obj = GlobalSkyModel(freq_unit,basemap,interpolation)
            % GLOBALSKYMODEL class constructor method.
            % obj = GlobalSkyModel(freq_unit,basemap,interpolation)
            %
            % Notes
            %     -----
            %     The matlab `spline` function allows one to explicitly
            %     set second derivatives to zero at the endpoints, as is done in
            %     the original GSM. We force this, but results might differ.
            %     Note that we use 'cubic' as input for compatibility with 
            %     pyGSM, but 'spline' is called internally. Further, we default
            %     to use PCHIP interpolation as in pyGSM, which does not set zero end slopes.
            % 
            % Inputs
            % - freq_unit: {'Hz',('MHz'),'GHz'}
            % - basemap: {'haslam','wmap',('5deg')}
            %             GSM version to generate. The 5deg map has 5.1 degree resolution.
            %             This is a synthesized map made of all the maps considered in
            %             the de Oliveira-Costa et. al. paper
            %             At frequencies below 1GHz, haslam is preferred; you get higher
            %             resolution (1 degree) by locking to the Haslam 408 MHz map.
            %             At CMB frequencies, it is best to lock to the WMAP 23 GHz
            %             map, which is presented denoised with 2 degree resolution.
            % - interpolation: {('pchip'),'cubic'}
            %             Choose whether to use cubic spline interpolation or
            %             piecewise cubic hermitian interpolating polynomial (PCHIP).
            %             PCHIP is designed to never locally overshoot data, whereas
            %             splines are designed to have smooth first and second derivatives.
            %
            % Outputs
            % - obj:    GlobalSkyModel object
            %
            % Dependencies
            % -
            %
            % Created: 2021-05-13, Dirk de Villiers
            % Updated: 2021-05-13, Dirk de Villiers
            %
            % Tested : Matlab R2020a
            %  Level : 1
            %   File : 
            %
            % Example
            %   GSM = GlobalSkyModel;
            
            warning('This class is a WIP: not all functionality done yet')
            
            if nargin > 0 && ~isempty(freq_unit), obj.freq_unit = freq_unit; end
            if nargin > 1 && ~isempty(basemap), obj.basemap = basemap; end
            if nargin > 2 && ~isempty(interpolation), obj.interpolation = interpolation; end
            
            % Set the path
            P = fileparts(mfilename('fullpath'));
            obj.dataPath = [P,'\data\gsm_components.h5'];
            
            obj = obj.update_interpolants;
        end
        
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
        
        
        function obj = update_interpolants(obj)
            
            % Load the PCA map from the .h5 file
            switch obj.basemap
                case '5deg'
                    compName = '/component_maps_5deg';
                case 'haslam'
                    compName = '/component_maps_408locked';
                case 'wmap'
                    compName = '/component_maps_23klocked';
            end
            obj.pca_map_data = transpose(h5read(obj.dataPath,compName));
            
            % Now, load the PCA eigenvalues
            pca_table = transpose(h5read(obj.dataPath,'/components'));
            obj.pca_freqs_mhz = pca_table(:,1);
            obj.pca_scaling   = pca_table(:,2);
            obj.pca_comps     = pca_table(:,3:end);
            
            % Interpolate to the desired frequency values
            ln_pca_freqs = log(obj.pca_freqs_mhz);
            
            switch obj.interpolation
                case 'cubic'
                    spl_scaling = spline(ln_pca_freqs, log(obj.pca_scaling));
                    spl1 = spline(ln_pca_freqs, [0;obj.pca_comps(:,1);0]);
                    spl2 = spline(ln_pca_freqs, [0;obj.pca_comps(:,2);0]);
                    spl3 = spline(ln_pca_freqs, [0;obj.pca_comps(:,3);0]);
                case 'pchip'
                    spl_scaling = pchip(ln_pca_freqs, log(obj.pca_scaling));
                    spl1 = pchip(ln_pca_freqs, obj.pca_comps(:,1));
                    spl2 = pchip(ln_pca_freqs, obj.pca_comps(:,2));
                    spl3 = pchip(ln_pca_freqs, obj.pca_comps(:,3));
            end
            obj.interp_comps = {spl_scaling,spl1,spl2,spl3};
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
            assert(min(freqs_mhz) >= 10 && min(freqs_mhz) <= 94000, 'Frequency values lie outside 10 MHz < f < 94 GHz')
            
            % Load interpolators and do interpolation
            ln_freqs = log(freqs_mhz);
            
            scaling = exp(ppval(obj.interp_comps{1},ln_freqs));
            spl1 = ppval(obj.interp_comps{2},ln_freqs);
            spl2 = ppval(obj.interp_comps{3},ln_freqs);
            spl3 = ppval(obj.interp_comps{4},ln_freqs);
            obj.generated_interp_comps = {scaling,spl1,spl2,spl3};
            
            spl = [spl1(:),spl2(:),spl3(:)];
            
            Nf = numel(freqs);
            map_out = zeros(size(obj.pca_map_data,1),Nf);
            for ff = 1:Nf
                map_out(:,ff) = sum(bsxfun(@times,obj.pca_map_data,spl(ff,:)),2).*scaling(ff);
            end
            obj.generated_map_data = map_out;
            obj.generated_map_freqs = freqs;
            
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
               gmap = obj.generated_map_data(:,idx);
               freq = obj.generated_map_freqs(idx);
           else
               gmap = obj.generated_map_data(:,1);
               freq = obj.generated_map_freqs(1);
           end
           
           if nargin < 3, logged = false; end
           if logged, gmap = log2(gmap); end
           
           healpixPlotMollweide(gmap)
           title(['Global Sky Model at ', num2str(freq), ' MHz from the ', obj.basemap, ' map'])
        end
        
        function plotCompsInterp(obj)
            
            obj = obj.generate(logspace(1,4.99,71));
            
            figure
            subplot 211
            semilogx(obj.pca_freqs_mhz,obj.pca_scaling,'ko'), grid on, hold on
            plot(obj.generated_map_freqs.*obj.freqScale,obj.generated_interp_comps{1},'k.-')
            ylabel('PCA scaling')
            title(['Interpolation: ',obj.interpolation])
            subplot 212
            color = 'kbr';
            for pp = 1:3
                semilogx(obj.pca_freqs_mhz,obj.pca_comps(:,pp),[color(pp),'o']), grid on, hold on
                plot(obj.generated_map_freqs.*obj.freqScale,obj.generated_interp_comps{pp+1},[color(pp),'.-'])
            end
            xlabel('Frequency (MHz)')
            ylabel('PCA components')
        end
    end
end
