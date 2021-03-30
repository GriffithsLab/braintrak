function import_raw_eeg
	% This function loads the mat file corresponding to the EDF data i.e.
	% e.g. psg_data\psgs\control_apnea\control_apnea_9\control_apnea_9.mat
	% It then performs Fourier transforms and artifact rejection, and then 
	% writes the output to psg_data/sleep_tfs
    
	idx = 1;

	for j = idx
		import('~/git/sleepmod/braintrak/edf_data/edf_1/edf_1.mat');
	end

function import(data_set)
	infile = fullfile(data_set);   %(os_prefix,sprintf('%s_%d',data_set,k),sprintf('%s_%d.mat',data_set,k));
	fprintf('Loading: %s\n',infile);
	fhandle = load(infile);
	colheaders = fhandle.colheaders;
    % fprintf('%s_%d: %s electrode\n',data_set,k,colheaders);
	[t,f,s,nspec,n_reject] = get_tfs(fhandle.t,fhandle.data,30,4,false);


	t_max_retain = floor(length(t)/30)*30;
	t = t(1:t_max_retain);
	t = t(:);
	s = s(:,1:t_max_retain,:);
	nspec = nspec(1:t_max_retain);

	state_str = cell(length(t),1);
	state_score = nan(1,length(t));

	states_available = {'W','5','1','2','3','4'};
	state_color_index = [2 3 4 5 6 6]; % Indexes for cdata

    state_str = cell(length(t),1);
    [state_str{:}] = deal('edf_x'); 
    state_score = fhandle.epoch_score;
    
% 	for j = 1:length(t)
% 		if mod(j-1,30) == 0 % This is on a 30s block
% 			state_str{j} = int2str(fhandle.epoch_score{(j-1)/30+1});
% 			dominant_state = fhandle.epoch_score{(j-1)/30+1};
%         else
% 			state_str{j} = sprintf('%s-%s (%d/%d)',fhandle.epoch_score{floor((j-1)/30)+1},fhandle.epoch_score{ceil((j-1)/30)+1},30-mod(j-1,30),mod(j-1,30));
% 			if 30-mod(j-1,30) > mod(j-1,30) % The dominant state is whichever one contributed more (time-wise)
% 				dominant_state = fhandle.epoch_score{floor((j-1)/30)+1};
% 			else
% 				dominant_state = fhandle.epoch_score{ceil((j-1)/30)+1};
% 			end
% 		end
% 		state_score(j) = state_color_index(strcmp(dominant_state,states_available));
% 	end
	outfile = fullfile('/home/taha/git/sleepmod/braintrak/edf_data/sleep_tfs');
	fprintf('Saving: %s\n',outfile);
	save(outfile,'t','f','s','state_score','state_str','nspec','n_reject','colheaders')
	fprintf('Done\n\n')
    
    
    
    
function [tv,fv,spectra,nspec,nreject] = get_tfs(t,V,window_length,fft_length,remove_artifacts)
		% Calculate t,f,s based on input time and voltage
		% nspec is the number of spectra that were averaged 
		if nargin < 5 || isempty(remove_artifacts)
			remove_artifacts = true;
		end

		if nargin < 4 || isempty(fft_length)
			fft_length = 4;
		end

		if nargin < 3 || isempty(window_length) % in seconds
			window_length = 30;
		end

		if max(abs(V)) > 500 % Some of the opioid data is not clipped to anything sensible
			clipping = 9999; % Disable the clipping check, since the data is just wierd
		else
			clipping = max(abs(V)) - 0.003; % (V is in mV) Clip if it is within 3uV of the maximum. 155 for control, 125 for some of the opioid
		end

		rate = round(1/(t(2)-t(1))); % The sampling frequency. window*rate is the number of points needed
		time_step = 1; % 1 second difference between successive spectra
		
		stoptime = floor(t(end))-fft_length; % (time in seconds) The time the last fft window starts
		[fv,~,P] = utils.rfft(V(1:rate*fft_length),rate);
		P = zeros(length(fv),stoptime+1);
		n_avg = (window_length-fft_length)/time_step+1; % Number of FFTs contributing to each spectrum

		clean = true(size(0:time_step:stoptime)); % Spectra all start out clean
		v_std = zeros(size(clean));
		% Get the set of 4s overlapping windows
		for j = 0:time_step:stoptime
			time_filter{j+1} = (1+(rate*j:(j+fft_length)*rate-1));
			[~,~,P(:,j+1)] = utils.rfft(V(time_filter{j+1}),rate);

			v_std(j+1) = std(V(time_filter{j+1}));

		end

		if remove_artifacts
			% Check clipping voltage
			clean_clipping = true(size(time_filter));
			for j = 1:length(time_filter)
				if sum(abs(V(time_filter{j}))>clipping) > 9 % Empirically, 9 is a good balance. This should be justified properly at some point
					clean_clipping(j) = false;
				end
			end

			% Check power standard deviations
			run_threshold = rate/2; % at least this many consecutive points for it to be an artifact
			in_range = @(pow,std_lim) pow < (mean(pow)+std(pow)*std_lim);
			hf_filter = fv > 30 & fv < 45; % hf cutoff is 4.5 Hz
			hf_pow = arrayfun(@(j) utils.mex_trapz(fv(hf_filter),P(hf_filter,j)),1:size(P,2));
			clean_range = in_range(v_std,4.5) & in_range(hf_pow,3); 

			% Check for runs of flat signals
			repeated_index = find(diff(V)==0); % The indexes that where i+1 = i
			runs = zeros(length(repeated_index),2);
			runs(1,:) = [repeated_index(1) repeated_index(1)];
			pointer = 1;
			for j = 2:length(repeated_index)
				if repeated_index(j) == runs(pointer,2)+1 % If this follows a run
					runs(pointer,2) = repeated_index(j);
				else
					pointer = pointer + 1;
					runs(pointer,:) = [repeated_index(j) repeated_index(j)];
				end
			end
			runs = runs(1:pointer,:);
			runs = [runs,runs(:,2)-runs(:,1)];
			runs = runs(runs(:,3)>run_threshold,:);
			run_reject = false(size(V));
			for j = 1:size(runs,1)
				run_reject(runs(j,1):runs(j,2)) = true;
			end

% 			figure
% 			plot(1:length(V),V)
% 			for j = 1:size(runs,1)
% 				set(gca,'XLim',[runs(j,1)-rate*3 runs(j,1)+rate*7])
% 				set(gca,'XTick',[runs(j,1)-rate*3:rate:runs(j,1)+rate*7],'XTick',[runs(j,1)-rate*3:rate:runs(j,1)+rate*7],'XTickLabel',[0:10])
% 				hold on
% 				asdf = plot(runs(j,1):runs(j,2),V(runs(j,1):runs(j,2)),'r')
% 				pause
% 				delete(asdf)
% 			end

			% Now go over the 4s windows again and see if they contained any run artifacts

			clean_runs = true(size(time_filter));
			for j = 1:length(time_filter)
				if any(run_reject(time_filter{j}))
					clean_runs(j) = false;
				end
			end
			
			clean = clean_clipping & clean_range & clean_runs;
			
		else
			clean = true(size(delta_pow));
		end

		% Now do the calculations for the 30 second spectra
		spectra = nan(length(fv),size(P,2)-n_avg);
		tv = 1:(size(P,2)-n_avg);
		
		nreject = sum(~clean);
		nspec = zeros(size(tv));

		for j = tv
			time_filter = j:(j+n_avg-1);
			time_filter = time_filter(clean(time_filter));
			nspec(j) = length(time_filter);
			if nspec(j) > 0
				spectra(:,j) = sum(P(:,time_filter),2);
			end
		end
		
		% And finally select only some of the frequencies and renormalize
		ffilt = fv > 0 & fv < 45;
		fv = fv(ffilt);
		spectra = spectra(ffilt,:);
		spectra = bsxfun(@rdivide,spectra,trapz(fv,spectra));

