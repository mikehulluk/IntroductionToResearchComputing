classdef Neurons
	
	properties
		dimensionality = uint32(1);
		preferredStimulus = 0;
		popSize = uint32(1);
		maxRate = 0;
		backgroundRate = 0;
		integrationTime = 0;
		distribution = 'Gaussian';
		a = 1.0;
		alpha = 0.5;
		R = [];
		add = 0.0;
		exponent = 1.0;
		truncate = false;
	end
	
	methods
		
		function obj = Neurons(varargin)
		%	NEURONS Neuron population class constructor.
		%	n = Neurons(dimensionality, preferredStimulus, maxRate, backgroundRate, integrationTime, variabilityScheme, variabilityOpts)
		%
		%	dimensionality - stimulus dimensionality (only 1-D stimuli currently supported)
		%	maxRate - maximum firing rate (Hz)
		%	backgroundRate - background (spontaneous) firing rate (Hz)
		%	variability - a Variability object
		%	integrationTime - spike counting time per trial
		%
		%	maxFiringRate and backgroundFiringRate can be scalars or vectors of length popSize

			%varargin
			%size(varargin{2})
                        
			switch nargin
			case 7
				% Standard constructor	
				if length(varargin{1}) == 1 && isnumeric(varargin{1})
					obj.dimensionality = uint32(varargin{1});
				else
					error([inputname(1) ' is not a valid stimulus dimensionality'])
				end
				
				if isnumeric(varargin{2}) & size(varargin{2}, 1) == obj.dimensionality
					obj.preferredStimulus = varargin{2}';
					obj.popSize = size(varargin{2}, 2);
				else
					error([inputname(2) ' is not a valid preferred stimulus value or vector'])
				end

				if length(varargin{3}) == 1 && isnumeric(varargin{3})
					obj.maxRate = double(varargin{3}(ones(obj.popSize, 1)));
				elseif length(varargin{3}) == obj.popSize && isvector(varargin{3}) && isnumeric(varargin{3})
					obj.maxRate = reshape(double(varargin{3}), obj.popSize, 1);
				else
					error([inputname(3) ' is not a valid maximum firing rate value or vector for population size ' obj.popSize])
				end

				if length(varargin{4}) == 1 && isnumeric(varargin{4})
					obj.backgroundRate = double(varargin{4}(ones(obj.popSize, 1)));
				elseif length(varargin{4}) == obj.popSize && isvector(varargin{4}) && isnumeric(varargin{4})
					obj.backgroundRate = reshape(double(varargin{4}), obj.popSize, 1);
				else
					error([inputname(4) ' is not a valid background firing rate value or vector for population size ' obj.popSize])
				end

				if length(varargin{5}) == 1 && isnumeric(varargin{5})
					obj.integrationTime = double(varargin{5});
				else
					error([inputname(5) ' is not a valid integration time'])
				end

				switch lower(varargin{6})
				case 'poisson'
					obj.distribution = 'Poisson';
					obj.a = [];
					obj.alpha = [];
					obj.R = [];
					obj.add = 0.0;
				case 'gaussian-independent'
					obj.distribution = 'Gaussian';
					obj.a = varargin{7}(1); % need checks
					obj.alpha = varargin{7}(2); % need checks
					obj.R = eye(obj.popSize);

					if length(varargin{7}) == 3
						obj.add = varargin{7}(3); % need checks
					else
						obj.add = 0.0;
					end
				case 'gaussian-uniform'
					obj.distribution = 'Gaussian';
					obj.a = varargin{7}(1); % need checks
					obj.alpha = varargin{7}(2); % need checks
					obj.R = varargin{7}(3) * ~eye(obj.popSize) + eye(obj.popSize);
					obj.add = 0.0;
				case 'gaussian-exponential'
					obj.distribution = 'Gaussian';
					obj.a = varargin{7}(1); % need checks
					obj.alpha = varargin{7}(2); % need checks
					c = varargin{7}(3); % need checks
					rho = varargin{7}(4); % need checks
					prefDiff = repmat(obj.preferredStimulus, 1, obj.popSize);
					prefDiff = prefDiff - prefDiff.';
					obj.R = c .* exp(-abs(double(prefDiff)) ./ rho) .* ~eye(obj.popSize) + eye(obj.popSize);
					obj.add = 0.0;
				case 'gaussian-gaussian'
					obj.distribution = 'Gaussian';
					obj.a = varargin{7}(1); % need checks
					obj.alpha = varargin{7}(2); % need checks
					c = varargin{7}(3); % need checks
					beta = 1.0 ./ degToRad(varargin{7}(4)).^2; % need checks
					prefDiff = repmat(obj.preferredStimulus, 1, obj.popSize);
					prefDiff = prefDiff - prefDiff.';
					obj.R = c .* exp((cosd(double(prefDiff)) - 1) .* beta) .* ~eye(obj.popSize) + eye(obj.popSize);
					obj.add = 0.0;
				case 'cercal'
					obj.distribution = 'Gaussian';
					obj.add = varargin{7}(1);
					obj.a = varargin{7}(2);
					obj.alpha = 0.5;
					obj.R = eye(obj.popSize);
					obj.exponent = 2.0;
				otherwise
					error([varargin{6} ' is not a valid variability regime'])
				end

			otherwise
				error('Wrong number of arguments')
			end
		end
		
		function i = mi(obj, method, stim, tol)
			if sum(strcmp(method, {'quadrature' 'randMC' 'quasirandMC'})) == 0
				error([method ' is not a valid SSI calculation method'])
			end

			if ~isa(stim, 'StimulusEnsemble')
				error([inputname(3) ' is not a SimulusEnsemble object'])
			end

			% obj.popSize x stim.n
			rMean = obj.integrationTime .* meanR(obj, stim);
			rMeanCell = squeeze(mat2cell(rMean, obj.popSize, ones(stim.n, 1)));

			% Compute mean response dependent cov matrix stack Q [ (popSize x popSize) x stim.n ]
			QCell1 = obj.Q(rMeanCell);

			% Compute Cholesky decomps of Q matrices
			cholQ = cellfun(@(q) chol(q)', QCell1, 'UniformOutput', false);
			
			% Invert Q matrices
			invQCell = cellfun(@inv, QCell1, 'UniformOutput', false);
			cholInvQCell = cellfun(@chol, invQCell, 'UniformOutput', false);

			% Replicate cell arrays
			cholInvQCell = repmat(cholInvQCell', [1 stim.n]);
			rMeanCella = repmat(rMeanCell', [1 stim.n]);
			
			% Define function for multivariate gaussian sampling
			% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
			
			if obj.truncate
                fRand = @(m, c, z) max((m + c * z), 0.0); % truncate
            else
                fRand = @(m, c, z) m + c * z; % don't truncate
            end
			
			iter = 1;
			miBuf = rand(100,1);
			varBuf = 0;
			cont = true;
			im = 0;
			MI = 0;

			while cont
				if ~mod(iter, 100)
					fprintf('MI  iter: %d  val: %.2g  var: %.4e\n', iter, MI, varBuf)
				end

				switch method
				case 'randMC'
					% Sample s from stimulus distribution
					[dummy bin] = histc(rand(), cumsum(stim.pS));
					bin = bin + 1;
					%s = double(stim.ensemble);
					%s = s(bin);

					% Sample r from response distribution
					% Generate vector of independent normal random numbers (mu=0, sigma=1)
					z = randn(obj.popSize, 1);
					% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
					% !!! NOTE NEGATIVE RESPONSES MAY BE TRUNCATED TO ZERO - see above !!!
					r = fRand(rMeanCell{bin}, cholQ{bin}, z);
				end

				% log P(r|s)
				% Replicate to form a stim.n x stim.n cell array of response vectors
				rCell = repmat({r}, [stim.n 1]);
				% Calculate response log probability densities
				lpRgS = cellfun(@normpdfln, rCell, rMeanCella(:,bin), cholInvQCell(:,bin), repmat({'inv'}, [stim.n 1]));

				% log P(r,s')
				% Mutiply P(r|s) and P(s) to find joint distribution
				pS = stim.pS;
				lpRS = lpRgS + log(pS');

				% log P(r)
				% Calculate marginal by summing over s'
				lpR = logsumexp(lpRS);
				
				% log P(s)
				lpS = log(pS(bin));

				% MI in bits (convert from log_e to log_2)
				lpRS = lpRS(bin);
				im = im + (lpRS - (lpR + lpS)) ./ log(2);

				% MI
				MI = im / iter;

				miBuf(1+mod(iter, 100)) = MI;
				varBuf = var(miBuf);
				cont = varBuf > tol & iter < 100000 | iter < 2000;

				iter = iter + 1;
			end

			i = MI;
		end
		
		function varargout = ssiss(obj, n, method, stim, stimOrds, tol, maxit, cont)
            try
                dummy = toc;
            catch
                tic
            end

			try
				% Test sanity of neuron indices
				obj.preferredStimulus(n);
			catch err
				error([inputname(2) ' is not a valid neuron index'])
			end
			
			try
				% Test sanity of stimulus ordinate indices
				stim.ensemble(stimOrds);
			catch err
				error([inputname(5) ' is not a valid stimulus ordinate index'])
			end
			
			if sum(strcmp(method, {'quadrature' 'randMC' 'quasirandMC'})) == 0
				error([method ' is not a valid SSI calculation method'])
			end

			if ~isa(stim, 'StimulusEnsemble')
				error([inputname(4) ' is not a SimulusEnsemble object'])
			end
			
			% Should we treat the function as continuous, i.e. can we use smoothness measures?
			switch lower(cont)
			case 'cont'
				continuous = true;
			case 'disc'
				continuous = false;
			otherwise
				error('Valid options are disc and cont')
			end
			
			% Create mask for calculating specific stimulus ordinates only
			if ~isempty(stimOrds)
				sMask = false(stim.n, 1);
				sMask(stimOrds) = true;
				sMaskN = sum(sMask + 0);
			else
				sMask = true(stim.n, 1);
				sMaskN = stim.n;
				stimOrds = 1:stim.n;
			end
			
			% Get mean responses for each stimulus ordinate
			% obj.popSize x stim.n
			rMean = obj.integrationTime .* meanR(obj, stim);
			rMeanCell = squeeze(mat2cell(rMean, obj.popSize, ones(stim.n, 1)));

			% Compute mean response dependent cov matrix stack Q [ (popSize x popSize) x stim.n ]
			QCell1 = obj.Q(rMeanCell);

			% Compute Cholesky decomps of Q matrices
			cholQ = cellfun(@(q) chol(q)', QCell1, 'UniformOutput', false);

			% Invert Q matrices and compute Cholesky decomps
			invQCell = cellfun(@inv, QCell1, 'UniformOutput', false);
			cholInvQCell = cellfun(@chol, invQCell, 'UniformOutput', false);
			
			if ~isempty(n)
				% Create logical vector (mask) identifying neurons that are *not* part of the marginal SSI
				margMask = ones(obj.popSize, 1);
				margMask(n) = false;
				% Number of remaining neurons
				margMask = logical(margMask);
				
				% Get mean responses for each stimulus ordinate
				rMeanMargCell = cellfun(@(r) r(margMask), rMeanCell, 'UniformOutput', false);
				
				% Compute mean response dependent cov matrix stack Q
				QCellMarg1 = cellfun(@(q) q(margMask, margMask), QCell1, 'UniformOutput', false);
								
				% Invert Q matrices and compute Cholesky decomps
				invQCellMarg = cellfun(@inv, QCellMarg1, 'UniformOutput', false);
				cholInvQCellMarg = cellfun(@chol, invQCellMarg, 'UniformOutput', false);
			end
			
			% Define function for multivariate gaussian sampling
			% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
			% Comment as appropriate if you want to truncate at zero
			% This will mess up the Gaussianity
			
			if obj.truncate
                fRand = @(m, c, z) max((m + c * z), 0.0); % truncate
            else
                fRand = @(m, c, z) m + c * z; % don't truncate
            end
			
			iter = 1;
			cont = true;
			hfPwrSSI = 0;
			iSP = zeros(1, sMaskN);
			iSPmarg = iSP;
			hfPwrIsur = 0;
			acc = zeros(1, sMaskN);
			accMarg = acc;
			
			while cont
				if ~mod(iter, 10)
					if continuous
						hfPwrSSI = hfpwr1(SSI);
						hfPwrIsur = hfpwr1(Isur);
						fprintf('SSISS iter: %d  HF power: %.4e %.4e\n', iter, hfPwrSSI, hfPwrIsur)
					else
						fprintf('SSISS iter: %d of %d\n', iter, maxit)
					end
				end

				switch method
				case 'randMC'
					% Sample r from response distribution
					% Generate vector of independent normal random numbers (mu=0, sigma=1)
					zCell = mat2cell(randn(obj.popSize, sMaskN), obj.popSize, ones(sMaskN, 1));
					% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
					% !!! NOTE NEGATIVE RESPONSES MAY BE TRUNCATED TO ZERO, SEE ABOVE !!!
					rCell = cellfun(fRand, rMeanCell(sMask), cholQ(sMask), zCell, 'UniformOutput', false); % stim.n cell array of obj.popSize vectors

				case 'quasirandMC'
					% 

				end

				% log P(r|s')
				% Calculate response probability densities
				lpRgS = cell2mat(cellsxfun(@normpdfln, rCell, rMeanCell', cholInvQCell', {'inv'}));
                
				% log P(r,s')
				% Mutiply P(r|s) and P(s) to find joint distribution
				lpRS = bsxfun(@plus, lpRgS, log(stim.pS')); % stim'.n x stim

				% log P(r)
				% Calculate marginal by summing over s'
				lpR = logsumexp(lpRS);
				
				% log P(s'|r)
				% Divide joint by marginal P(r)
				lpSgR = bsxfun(@minus, lpRS, lpR);

				% H(s'|r), in bits, converting from log_e to log_2
				hSgR = -sum(exp(lpSgR) .* (lpSgR ./ log(2)), 1);

				% Accumulate Isp(r)
				% Specific information; reduction in stimulus entropy due to observation of r
				iSP = iSP + stim.entropy - hSgR;
				
				% log_2( P(r|s) / P(r) )
				% Accumulate samples
				acc = acc + (diag(lpRgS(stimOrds,:))' - lpR) ./ log(2);
				
				% Issi(s)
				% SSI; average of Isp over all r samples
				fullSSI = iSP ./ iter;
				
				% Isur(s)
				fullIsur = acc ./ iter;
				
				if exist('margMask', 'var')
					% If we are calculating a marginal SSI, compute the SSI for remaining neurons

					% Mask out neurons of interest in response vectors
					rCellMarg = cellfun(@(r) r(margMask), rCell, 'UniformOutput', false);

					% log P(r|s)
					lpRgS = cell2mat(cellsxfun(@normpdfln, rCellMarg, rMeanMargCell', cholInvQCellMarg', {'inv'}));

					% log P(r,s)
					% Multiply P(r|s) and P(s) to find joint distribution
					lpRS = bsxfun(@plus, lpRgS, log(stim.pS')); % stim'.n x stim

					% log P(r)
					% Calculate marginal by summing over s'
					lpR = logsumexp(lpRS);

					% log P(s|r)
					% Divide joint by marginal P(r)
					lpSgR = bsxfun(@minus, lpRS, lpR);

					% H(s|r), in bits, converting from log_e to log_2
					hSgR = -sum(exp(lpSgR) .* (lpSgR ./ log(2)), 1);

					% Isp(r)
					% Specific information; reduction in stimulus entropy due to observation of r
					iSPmarg = iSPmarg + stim.entropy - hSgR;
					
					% log_2( P(r|s) / P(r) )
					% Accumulate samples
					accMarg = accMarg + (diag(lpRgS(stimOrds,:))' - lpR) ./ log(2);

					% Issi(s)
					% SSI; average of Isp over all r samples
					remainderSSI = iSPmarg ./ iter;
					SSI = fullSSI - remainderSSI;
					
					% Isur(s)
					remainderIsur = accMarg ./ iter;
					Isur = fullIsur - remainderIsur;
				else
					remainderSSI = fullSSI;
					SSI = fullSSI;
					
					remainderIsur = fullIsur;
					Isur = fullIsur;
				end

				% Smoothness measure doesn't work if we are only
				% calculating selected stimulus ordinates, so run until
				% iteration limit
				if ~mod(iter, 10) && continuous
					cont = hfPwrSSI > tol || hfPwrIsur > tol;
				else
					cont = true;
				end

				if iter < 100
					cont = true;
				end

				if iter >= maxit
					cont = false;
					disp('Iteration limit exceeded')
				end
				
				if exist('/tmp/haltnow', 'file')
					cont = false;
					disp('Detected /tmp/haltnow, aborting calculation')
				end
				
				% Halt before eddie kills the job
				if (toc / 3600.0) > 47.75
					cont = false;
					disp('Runtime approaching 48 hrs, halting calculation')
				end

				iter = iter + 1;
			end

			fprintf('SSISS iter: %d  elapsed time: %.4f seconds\n', iter - 1, toc)

			switch nargout
			case 1
				varargout = {SSI};
			case 2
				varargout = {SSI (iter - 1)};
			case 3
				varargout = {fullSSI remainderSSI (iter - 1)};
			case 5
				varargout = {fullSSI remainderSSI fullIsur remainderIsur (iter - 1)};
			otherwise
				error('Wrong number of outputs')
			end
		end
		
		function varargout = fisher(obj, method, stim, tol, maxit)
			
			if ~isa(stim, 'StimulusEnsemble')
				error([inputname(4) ' is not a SimulusEnsemble object'])
            end
            
            try
                dummy = toc;
            catch
                tic
            end
            
			switch method
			case 'analytic'	
				% Info=Fisher(N, alpha, correlation_type, param_correlation)
				% compute I_mean, I_cov and I_mean and I_cov if neurons are independent 
				% when the tuning function is circular normal and correlation is toeplitz
				% 
				% PARAMETERS:
				% -----------
				% - N is the number of neurons
				% - f is the tuning curve
				% - f_prime is the derivative of the tuning curve
				% - Q is the covariance matrix
				% - k is the fano factor, now a vector
				% - k prime is its derivative
				%
				% RETURNS:
				% ---------
				% Info(1):  I_mean
				% Info(2):  I_mean_ind if the correlations were set to 0 (same variance)
				% Info(3):  I_cov
				% Info(4):  I_cov_ind if the correlations were set to 0
				% Info(5):  I_mean+I_cov
				% Info(6):  I_SQE

				% pseries@gatsby.ucl.ac.uk  2/20/2005

				%phi = repmat(degToRad(double(obj.preferredStimulus)) [1 stim.n]);
				%phi_stim = repmat(degToRad(double(stim.ensemble)) [obj.popSize 1]);
				%x = phi - phi_stim; % in radians, phi_stim is the orientation of the stimulus

				% ==============================================================
				% Tuning Curves (circular normal), mean variance relationship, and
				% useful Fourier transforms
				% ==============================================================

				f = obj.integrationTime .* meanR(obj, stim);
				f_prime = obj.integrationTime .* dMeanR(obj, stim);

				g0 = f_prime ./ f;
				%h0 = f_prime ./ (f.^obj.alpha);

				%g0_tilda= ifftshift(abs(ifft(g0)));  % and Fourier
				%tranforms
				%h0_tilda_square= conj(ifftshift(ifft(h0))).*ifftshift(ifft(h0));
				%g0_tilda_square= conj(ifftshift(ifft(g0))).*ifftshift(ifft(g0));

				% ==============================================================
				% Correlation matrix, Covariance matrix, inverse and derivative.
				% Fourier tranform (eigenvalues)
				% ==============================================================

				fCell = squeeze(mat2cell(f, obj.popSize, ones(stim.n, 1)));
				f_primeCell = squeeze(mat2cell(f_prime, obj.popSize, ones(stim.n, 1)));
				g0Cell = squeeze(mat2cell(g0, obj.popSize, ones(stim.n, 1)));
				QCell1 = obj.Q(fCell);
				Q_inv = cellfun(@inv, QCell1, 'UniformOutput', false); % inverse
				k = repmat({obj.a(ones(obj.popSize,1))}, [1 stim.n]);
				k_prime = repmat({zeros(obj.popSize,1)}, [1, stim.n]);
				
				fQ_prime = @(q, kk, kp, gz) obj.alpha * (diag(gz) * q + q * diag(gz)) + (diag(kp ./ kk)) * q;
				Q_prime = cellfun(fQ_prime, QCell1, k, k_prime, g0Cell, 'UniformOutput', false); % derivative

				%d = cellfun(@(kk, ff) kk .* ff.^(2*obj.alpha), k, fCell, 'UniformOutput', false);	
				%d_prime = cellfun(@(kk, kp, ff, fp) 2 .* obj.alpha .* kk .* fp .* ff.^(2*obj.alpha-1) + kp .* ff.^(2*obj.alpha), k, k_prime, fCell, f_primeCell, 'UniformOutput', false);
				%D_inv = cellfun(@(dd) inv(diag(dd)), d, 'UniformOutput', false);
				%D_prime = cellfun(@diag, d_prime, 'UniformOutput', false);
				%J = cellfun(@times, QCell1, QCell1, 'UniformOutput', false);

				% ==============================================================
				% Fisher Information
				% I_mean, I_cov and approximation of I_cov
				% ==============================================================

				if obj.popSize > 1
					Info(1,:) = cellfun(@(fp, qi) fp' * qi * fp, f_primeCell, Q_inv);         % direct method
					%Info(2) = cellfun(@(fp, di) fp * di * fp', f_primeCell, D_inv);
				else
					Info(1,:) = cellfun(@(fp, q) fp.^2 / q, f_primeCell, QCell1);
				end

				Info(3,:) = cellfun(@(qp, qi) 0.5 * trace(qp * qi * qp * qi), Q_prime, Q_inv);  % direct method for Icov
				%Info(4) = 0.5 * trace(D_prime * D_inv * D_prime * D_inv);  % direct method for Icov_shuffle

				Info(5,:) = Info(1,:) + Info(3,:);                         % Fisher
				%Info(6) = 0.5 * d_prime * inv(J) * d_prime';             % second term of I_sqe

				%Info = Info * (pi^2) / (180.^2); % convert to degrees ^-2
				%fprintf('\n.. I_mean : %.5f \n.. I_mean_ind : %.5f \n.. I_cov : %.5f \n.. I_cov_ind : %.5f \n.. I_tot : %.5f\n\n', Info(1), Info(2), Info(3), Info(4), Info(5));
				%fprintf('I_SQE : %.5f\n\n', Info(1)+Info(6));

				switch nargout
				case 1
					varargout = {Info(5,:)};
				case 2
					varargout = {Info(1,:) Info(3,:)};
				otherwise
					error('Wrong number of outputs')
				end

			case {'randMC' 'quasirandMC'}
				% obj.popSize x stim.n
				rMean = obj.integrationTime .* meanR(obj, stim);
                rMeanCell = squeeze(mat2cell(rMean, obj.popSize, ones(stim.n, 1)));

				% Compute mean response dependent cov matrix stack Q [ (popSize x popSize) x stim.n ]
                QCell1 = obj.Q(rMeanCell);
                
                % Compute Cholesky decomps of Q matrices
                cholQ = cellfun(@(q) chol(q)', QCell1, 'UniformOutput', false);
                
                % Invert Q matrices and compute Cholesky decomps
                invQCell = cellfun(@inv, QCell1, 'UniformOutput', false);
                cholInvQCell = cellfun(@chol, invQCell, 'UniformOutput', false);
                
				% Normalisation terms for multivariate Gaussian
				%normFactor = 1.0 ./ ((2.0 .* pi).^(double(obj.popSize) ./ 2.0) .* cellfun(@det, QCell1).^0.5);
				%normFactor = mat2cell(normFactor, 1, ones(length(normFactor), 1));
				%normFactor = num2cell(normFactor);

				% Replicate normalisation factors and cov matrices into a 3-row cell array.  This allows us to
				% compute two d/ds values (one either side of the nominal s) for each s.
				%normFactor = [normFactor(end) normFactor(1:end-1) ; normFactor ; normFactor(2:end) normFactor(1)];
				QCell = [QCell1(end) QCell1(1:end-1) ; QCell1 ; QCell1(2:end) QCell1(1)];

				% Define function for multivariate gaussian sampling
				% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
                
				if obj.truncate
				    fRand = @(m, c, z) max((m + c * z), 0.0); % truncate
				else
				    fRand = @(m, c, z) m + c * z; % don't truncate
                end
                
				iter = 1;
                cont = true;
                hfPwrFI = 0;
                acc = zeros(1, stim.n);
                
                while cont
                    if ~mod(iter, 10)
                        hfPwrFI = hfpwr1(FI);
                        fprintf('Fisher iter: %d  HF power: %.4e\n', iter, hfPwrFI)
                    end
                    
					switch method
					case 'randMC'
						% Sample r from response distribution
						% Generate vector of independent normal random numbers (mu=0, sigma=1)
						zCell = squeeze(mat2cell(randn(obj.popSize, stim.n), obj.popSize, ones(stim.n, 1)));
						% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
						% !!! NOTE NEGATIVE RESPONSES ARE TRUNCATED TO ZERO !!!
						rCell = cellfun(fRand, rMeanCell, cholQ, zCell, 'UniformOutput', false); % stim.n cell array of {obj.popSize}
                        
					case 'quasirandMC'
						%
                    end

                    % log P(r|s)
                    % Calculate response probability densities
                    lpRgS = cell2mat(cellsxfun(@normpdfln, rCell, rMeanCell, cholInvQCell, {'inv'})); % 1 x stim.n
                    
                    % d/ds log P(r|s)
					% Find numerical gradient
					dlpRgS = gradient(lpRgS); % 1 x stim.n
                    
					% Divide by bin widths to get derivative
					dlpRgS = dlpRgS ./ stim.width; % 1 x stim.n

					% (d/ds log P(r|s))^2
					acc = acc + dlpRgS .^ 2;

					switch method
					case 'randMC'
						%FI = mean(dlpRgS2, 1);
						FI = acc ./ iter;
					end

                    % Smoothness measure doesn't work if we are only
                    % calculating selected stimulus ordinates, so run until
                    % iteration limit
                    if ~mod(iter, 10)
                        cont = hfPwrFI > tol;
                    else
                        cont = true;
                    end
    
                    if iter < 100
                        cont = true;
                    end
    
                    if iter >= maxit
                        cont = false;
                        disp('Iteration limit exceeded')
                    end
                    
                    if exist('/tmp/haltnow', 'file')
                        cont = false;
                        disp('Detected /tmp/haltnow, aborting calculation')
                    end
                    
                    % Halt before eddie kills the job
                    if (toc / 3600.0) > 47.75
                        cont = false;
                        disp('Runtime approaching 48 hrs, halting calculation')
                    end
    
                    iter = iter + 1;
				end

				varargout = {FI};
			otherwise
				error([method ' is not a valid FI calculation method'])
			end
		end
		
		function varargout = margFisher(obj, nMarg, method, stim, tol, varargin)
			if isempty(varargin)
				option = 'raw';
			else
				option = varargin{1};
			end

			if ~isa(stim, 'StimulusEnsemble')
				error([inputname(4) ' is not a SimulusEnsemble object'])
			end

			% Compute FI of full population
			fullFisher = fisher(obj, method, stim, tol);
			
			% Create population of remaining neurons
			remainingNeurons = obj.remove(nMarg);
			
			% FI of remaining neurons
			remainderFisher = fisher(remainingNeurons, method, stim, tol);
			
			switch option
			case 'diff'
				varargout = {fullFisher - remainderFisher};
			case 'rootDiff'
				varargout = {sqrt(fullFisher - remainderFisher)};
			case 'diffRoot'
				varargout = {sqrt(fullFisher) - sqrt(remainderFisher)};
			case 'invDiffInv'
				varargout = {fullFisher .* remainderFisher ./ (remainderFisher - fullFisher)};
			case 'diffSD'
				varargout = {(fullFisher.^-0.5 - remainderFisher.^-0.5).^-2};
			case 'raw'
				varargout = {fullFisher remainderFisher};
			otherwise
				error('Invalid marginal fisher option')
			end
		end
		
		function ifish = Ifisher(obj, stim)
			fish = obj.fisher('analytic', stim, 0);
			pS = stim.pS;
			zOrds = find(fish == 0);
			
			if ~isempty(zOrds)
				fish(zOrds) = [];
				pS(zOrds) = [];
				pS = pS ./ sum(pS);
			end
			
			ifish = stim.entropy - sum(pS .* 0.5 .* log2(2.0 .* pi .* exp(1) ./ fish));
		end
		
		function ifish = mIfisher(obj, nMarg, stim)
			FI = margFisher(obj, nMarg, 'analytic', stim, 0, 'raw');
			FItotal = FI(1,:);
			FIrem = FI(2,:);

			ifishTotal = stim.entropy - sum(stim.pS .* 0.5 .* log2(2.0 .* pi .* exp(1) ./ FItotal));
			ifishRem = stim.entropy - sum(stim.pS .* 0.5 .* log2(2.0 .* pi .* exp(1) ./ FIrem));

			ifish = ifishTotal - ifishRem;
		end
		
		function ssif = SSIfisher(obj, n, fisherMethod, stim, tol)
		
			% S stimulus
			sMat = repmat(stim.ensemble, [stim.n 1]);
			% sHat stimulus estimate
			sHatMat = repmat(stim.ensemble', [1 stim.n]);
			% log p(S)
			psMat = repmat(stim.pS, [stim.n 1]);
			
			if stim.circular
				% circular difference
				dS = mod(sHatMat - sMat, stim.circular);
				i = find(dS > (0.5 * stim.circular));
				dS(i) = dS(i) - stim.circular;
				i = find(dS < (-0.5 * stim.circular));
				dS(i) = dS(i) + stim.circular;
			else
				% linear difference
				dS = sHatMat - sMat;
			end
			
			if ~isempty(n)
				[fullFI remFI] = obj.margFisher(n, fisherMethod, stim, tol, 'raw');
				
				% Compute SSIfisher excluding cells of interest
				% sigma(s)
				% Compute SD of optimal estimator as a function of the stimulus
				sigma = remFI .^ -0.5;
				sigmaMat = repmat(sigma, [stim.n 1]);

				% log p(sHat|S)
				lpsHat_s = cellfun(@normpdfln, num2cell(dS), repmat({0}, [stim.n stim.n]), num2cell(sigmaMat));
				
				% log p(S)
				lpS = log(psMat);
				% log p(sHat,S)
				lpsHats = lpsHat_s + lpS;
				% log p(sHat)
				lpsHat = logsumexp(lpsHats, 2);
				% log p(S|sHat)
				lps_sHat = lpsHats - repmat(lpsHat, [1 stim.n]);

				% H(s|sHat) as a function of sHat
				Hs_sHat = -sum(exp(lps_sHat) .* lps_sHat ./ log(2), 2);

				% Isp(sHat) specific information
				isp = stim.entropy - Hs_sHat;

				% SSIfisher
				ssifRem = sum(exp(lpsHat_s) .* repmat(isp, [1 stim.n]), 1);
			else
				fullFI = obj.fisher(fisherMethod, stim, tol);
			end
			
			% Compute SSIfisher for full population
			% sigma(s)
			% Compute SD of optimal estimator as a function of the stimulus
			sigma = fullFI .^ -0.5;
			sigmaMat = repmat(sigma, [stim.n 1]);

			% log p(sHat|S)
			lpsHat_s = cellfun(@normpdfln, num2cell(dS), repmat({0}, [stim.n stim.n]), num2cell(sigmaMat));
			
			% log p(S)
			lpS = log(psMat);
			% log p(sHat,S)
			lpsHats = lpsHat_s + lpS;
			% log p(sHat)
			lpsHat = logsumexp(lpsHats, 2);
			% log p(S|sHat)
			lps_sHat = lpsHats - repmat(lpsHat, [1 stim.n]);

			% H(s|sHat) as a function of sHat
			Hs_sHat = -sum(exp(lps_sHat) .* lps_sHat ./ log(2), 2);

			% Isp(sHat) specific information
			isp = stim.entropy - Hs_sHat;
			
			% SSIfisher
			ssifFull = sum(exp(lpsHat_s) .* repmat(isp, [1 stim.n]), 1);
			
			if ~isempty(n)
				ssif = ssifFull - ssifRem;
			else
				ssif = ssifFull;
			end
		end
		
		function obj = gainadapt(obj, width, amnt, centre)
			obj.maxRate = obj.maxRate .* (1 - amnt .* exp(-(1.0 ./ degToRad(width).^2) .* (1 - cosd(double(obj.preferredStimulus - centre)))));
		end
		
		function p = pOfR(obj, response, stimulusSet)
			if ~(isnumeric(response) && isvector(response) && length(response) == obj.popSize)
				error([inputname(2) ' is not a valid response'])
			end

			if ~(isnumeric(stimulusSet) && size(stimulusSet, 1) == obj.dimensionality)
				error([inputname(3) ' is not a valid stimulus set'])
			end

			r = repmat(reshape(response, length(response), 1), 1, length(stimulusSet));
			rMean = meanR(obj, stimulusSet);
			r = r - rMean;
			
			prob = zeros(1,length(stimulusSet));
			
			for s = 1 : length(stimulusSet)
				Q = obj.varA * obj.varR .* (rMean(:,s) * rMean(:,s)') .^ obj.varAlpha;
				normFactor = 1.0 / ((2.0 * pi)^(obj.popSize / 2.0) * det(Q)^0.5);
				prob(1,s) = normFactor * exp(r(:,s)' * (Q \ r(:,s)));
			end

			p = prob;
		end
		
		function q = Q(obj, resp)
			%if obj.add == 0
			%	q = cellfun(@(r) obj.a .* obj.R .* (r * r').^obj.alpha, resp, 'UniformOutput', false);
			%else
			%	q = cellfun(@(r) (obj.add + obj.a .* obj.R .* (r * r').^obj.alpha).^2, resp, 'UniformOutput', false);
			%end
			
			q = cellfun(@(r) (obj.add .* obj.R + obj.a .* obj.R .* (r * r').^obj.alpha).^obj.exponent, resp, 'UniformOutput', false);
		end
		
		function retStr = char(obj)
			%	CHAR Text representation of Neurons object
			%	aString = char(object)

			if ~isscalar(obj)
				retStr = 'Array of Neurons objects';
				return
			end

			prefStr = char(obj.preferredStimulus);
			maxRateStr = char(obj.maxRate);
			backgroundRateStr = char(obj.backgroundRate);
			rStr = display(obj.R);

			formatStr = ['Population size: %.0f\n'...
						 'Stimulus dimensionality: %.0f\n'...
						 'Preferred stimuli:\n'...
						  prefStr '\n'...
						 'Maximum firing rate (Hz):\n'...
						  maxRateStr '\n'...
						 'Background firing rate (Hz):\n'...
						  backgroundRateStr '\n'...
						 'Integration time: %.3f s\n'...
						 'Variability distribution: %s\n'...
						 'Variability coefficient a: %.3f\n'...
						 'Variability exponent alpha: %.3f\n'...
						 'Correlation matrix:\n'];

			retStr = strvcat(sprintf(formatStr, obj.popSize, double(obj.dimensionality), obj.integrationTime, obj.distribution, obj.a, obj.alpha), rStr);
		end
		
		function retVal = tcplot(obj, stim)
			r = meanR(obj, stim);
			[stims ind] = sort(double(stim.ensemble));
			retVal = plot(stims, r(:,ind));
		end
		
		function varargout = remove(obj, nMarg)
			if ~isempty(nMarg)
				% Create logical vector (mask) identifying neurons that are *not* part of the marginal SSI
				margMask = ones(obj.popSize, 1);
				margMask(nMarg) = false;
				% Number of remaining neurons
				nMarg = sum(margMask);
				margMask = logical(margMask);
			else
				error('Must specify index of a neuron or neurons')
			end
			
			obj.preferredStimulus = obj.preferredStimulus(margMask);
			obj.popSize = uint32(nMarg);
			
			if length(obj.maxRate) > 1
				obj.maxRate = obj.maxRate(margMask);
			end
			
			if length(obj.backgroundRate) > 1
				obj.backgroundRate = obj.backgroundRate(margMask);
			end
			
			Rmask = logical((margMask+0) * (margMask+0)');
			obj.R = reshape(obj.R(Rmask), [nMarg nMarg]);
			
			switch nargout
			case 1
				varargout = {obj};
			case 2
				varargout = {obj margMask};
			otherwise
				error('Wrong number of outputs')
			end
					
		end
		
	end

end
