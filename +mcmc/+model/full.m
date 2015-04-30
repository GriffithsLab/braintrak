classdef full < mcmc.model.template
	% This model is set up to fit t0 keeping the prior fixed if the alpha peak is not present
	% Note that prepare_for_fit must be used before spectrum() can be used
	% Uniform priors and DB initial fit

	properties
		kmax
		k0
		k2u
		k2_volconduct

		stab_gamma_prefactor
		spec_gamma_prefactor

		stab_w
		spec_w
		
		p
		weights
		emg
		normalization_target
		stab_Mtot
		spec_Mtot
		compute_Mtot
	end


	methods
		function self = full() % Constructor
			self.name = 'model_full';
			self.n_params = 9;
			self.param_names = {'Gee','Gei','Gese','Gesre','Gsrs','Alpha','Beta','t0','EMGa'};
			self.param_symbols = {'G_{ee}','G_{ei}','G_{ese}','G_{esre}','G_{srs}','\alpha','\beta','t_0','A_{EMG}'};
			self.param_units = {'','','','','','s^{-1}','s^{-1}','ms',''};

			self.n_fitted = 9;
			self.skip_fit = zeros(1,self.n_fitted);

			self.initial_step_size = [0.4  0.4  1  1  0.2  5  40  0.005  0.001];
		   	self.limits = [ eps  -40        eps     -40      -14      10      100    0.075   0  ;...
		      		        20   -eps       40      -eps      -eps      100    800    0.14   1  ];
			
			self.kmax = 4;
			self.k0 = 10; % Volume conduction parameter
			self.p = model.params;
			self.p.phin = 1e-5;
			% Limits based on stability_3d
			% Somewhat tighter ranges
			% EMG ranges based on spectfit output generated by sacha, speciically spec10000113.fit
		end

		function p = p_from_params(self,fitted_params) % Map parameters to point
			p = model.params;
			p.gabcd = fitted_params(1:5);
			p.alpha(:) = fitted_params(6);
			p.beta(:) = fitted_params(7);
			p.t0 = fitted_params(8);
			p.taues = p.t0/2;
			p.tause = p.t0/2;
			p.emg_a = fitted_params(9);
		end
		
		function params = params_from_p(self,p) % Map point to parameters
			params = [p.gabcd p.alpha(1) p.beta(1) p.t0 p.emg_a];
		end

		function valid = validate_params(self,pars)
			valid = ~(pars(1)/pars(2) > -0.5 || (pars(1) + pars(2)) > 1 || pars(7)/pars(6) > 20 || any(pars > self.limits(2,:) | pars < self.limits(1,:) ));
		end

		function [chisq,P] = objective(self,pars) % Calculate the objective
			P = [];
			chisq = NaN;

		    if self.compute_Mtot % If this flag has been set to false, then these fields are already present
		        % Should be fine to let that just throw an undefined variable error later if this isn't the case
		       self.stab_Mtot = exp(1i*self.stab_w*pars(8));
		       self.spec_Mtot = exp(1i*self.spec_w*pars(8));
		    end

			stab_L = 1./((1-1i*self.stab_w/pars(6)).*(1-1i*self.stab_w/pars(7)));
		    spec_L = 1./((1-1i*self.spec_w/pars(6)).*(1-1i*self.spec_w/pars(7)));

		    stab_gamma_prefactor = self.stab_gamma_prefactor;
		    d=(stab_gamma_prefactor.*(1-stab_L.*pars(2))-stab_L.*pars(1)).*(1-stab_L.*stab_L.*pars(5))-stab_L.*stab_L.*self.stab_Mtot.*(pars(3) + stab_L.*pars(4));      

		    stab=d(1)>0 && ~any(real(d(2:end))<0 & imag(d(2:end)).*imag(d(1:end-1))<0);

		    if ~stab
		    	return
		    end

		    % And the spectrum
		    Jei_oneminus = 1-spec_L.*pars(2);
		    Jsrs_oneminus = 1-spec_L.*spec_L.*pars(5);

		    re2 = self.p.re.^2;

		    q2re2 = (self.spec_gamma_prefactor - 1./(Jei_oneminus).*(spec_L.*pars(1) + ((spec_L.*spec_L.*pars(3) + spec_L.*spec_L.*spec_L.*pars(4)).*self.spec_Mtot)./(Jsrs_oneminus)));
		    T_prefactor = spec_L.*spec_L.*self.p.phin./(Jei_oneminus.*Jsrs_oneminus);
		    % T prefactor doesn't have a exp(i*omega*t0/2) term because the modulus of this is 1
		    P = zeros(size(self.spec_w));

		    k2u = self.k2u;
		    k2_volconduct = self.k2_volconduct;

		    for j = 1:size(k2u,1)
		        P = P + k2u(j,2).*abs(T_prefactor./(k2u(j,1)*re2+q2re2)).^2 * k2_volconduct(j); % For use with the efficient way
		    end

		    P = P + 1e-12*pars(9)*self.emg; % Add the EMG component
		    
		    P = P./utils.mex_trapz(self.target_f(self.weights>0),P(self.weights>0));
		    P = self.normalization_target*P;
		    sqdiff = (abs(P-self.target_P)./self.target_P).^2; % This is the squared fractional difference
		    chisq = sum(sqdiff(:).*self.weights(:));
		end
		
		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			[f1,P1,idx,db_data] = db_fit.quick_fit(self,target_f,target_P);
			%loglog(target_f,target_P,f1,P1)
			a =  db_data.gab(idx,:);
			p = model.params(db_data.iswake(idx));

			initial_values =  [a(1) a(2) a(3)*a(4) a(3)*a(5)*a(7) a(5)*a(8) p.alpha(1) p.beta(1) p.t0 0];

			prior_pp = self.uniform_priors();
		end
		
		function set_cache(self,initial_values) % Calculate the cache - takes in skip_fit

			% Cache for the mcmc fit code
			% Lower resolution, dynamic alpha and beta
			Lx = 0.5; % linear dimensions of cortex (assume square)
			dk = 2*pi/Lx;
			m_rows = -self.kmax:self.kmax; 
			n_cols = -self.kmax:self.kmax; 
			[kxa,kya] = meshgrid(dk*m_rows,dk*n_cols); 
			k2 = kxa.^2+kya.^2;
			k2u = unique(k2(:));
			self.k2u = [k2u histc(k2(:),k2u)]; % First col is the k2 value, second is how many times it appeared
			self.stab_w=2*pi*linspace(0,100,5000); % High res array, 0 to 100Hz with 10000 points
			self.spec_w = self.target_f*2*pi;

			self.stab_w = self.stab_w(:);
			self.spec_w = self.spec_w(:);

			self.stab_gamma_prefactor = (1-1i*self.stab_w/self.p.gammae).^2;
			self.spec_gamma_prefactor = (1-1i*self.spec_w/self.p.gammae).^2;

			self.weights = self.get_weights(self.target_f);

			self.k2_volconduct = exp(-k2u(:,1)/self.k0^2);

			emg_f = 40;
			self.emg = (self.spec_w/(2*pi*emg_f)).^2./(1+(self.spec_w/(2*pi*emg_f)).^2).^2;
			
			if ~isempty(self.target_P)
				% If this cache file is being generated for a specific fit instance
				self.normalization_target = utils.mex_trapz(self.target_f(self.weights>0),self.target_P(self.weights>0));
			end


			if self.skip_fit(8) == 2 && ~isempty(initial_values) % If t0 is not being fitted, then the value is of course fixed
				self.stab_Mtot = exp(1i*self.stab_w*initial_values(8));
				self.spec_Mtot = exp(1i*self.spec_w*initial_values(8));
				self.compute_Mtot = false;
			else
				self.compute_Mtot = true; 
			end
		end
		
		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz(:,1) = params(:,1)./(1-params(:,2));
			xyz(:,2) = (params(:,3)+params(:,4))./(1-params(:,5))./(1-params(:,2));
			xyz(:,3) = -params(:,5).*params(:,6).*params(:,7)./((params(:,6)+params(:,7)).^2);
		end
		
		function prepare_for_fit(self,varargin)
			% argument list: target_f,target_P,initial_values,prior_pp,skip_fit
			if nargin < 6 || isempty(varargin{5})
				self.default_t0_skip_fit(varargin{1},varargin{2},8);
			else
				assert(length(varargin{5}) == self.n_params);
				self.skip_fit = varargin{5};
			end
			
			prepare_for_fit@mcmc.model.template(self,varargin{1:4});
		end

	end
end


