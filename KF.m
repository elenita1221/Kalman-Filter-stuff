function [xf_means,xa_means,Pxxas,vs] =...
  KF(KF_choice,KF_model_func,KF_exp_func,KF_meas_func,KF_ic_state_func,...
  KF_ic_params_func,KF_meas_addl_info, KF_ic_state_addl_info,KF_ic_params_addl_info,...
  KF_noise_addl_info,KF_time_info,KF_options_cell,problem_location)
  
  %  function [xps,xas,Pxxps,Pxxas,par_means,par_vars] =...
  %   KF(KF_choice,KF_model_func,KF_exp_func,KF_meas_func,KF_ic_params_func,...
  %   KF_ic_state_func,KF_noise_func,KF_ic_state_addl_info,...
  %   KF_ic_paras_addl_info,KF_noise_addl_info,KF_time_info,KF_options_cell)
  %  
  %  See Call_KF or Call_Cell_Migration_KF to see some examples of the
  %  program call.
  %   
  %  Inputs:
  %
  %  KF_choice-string option:  'En','SC','KLEn','KLSC','KF'
  %    more options may be added later, but for now these are all the
  %    options.
  %
  %  KF_model_func-function handle.  The function should take in three
  %    arguments, a state vector, a time step, and a parameter structure
  %    (consisting of parameters that are not allowed to be free/are set
  %    during the simulations) and return a prediction for the state vector
  %    at the next time step.
  %
  %  KF_exp_func-function handle.  The function should take in one
  %    argument, a structure with the time(s) and, if a spatial problem, 
  %    locations at which you want the experimental measurements and return
  %    the corresponding experimental measurements.
  % 
  %  KF_meas_func-function handle.  Takes state and turns it into
  %    the corresponding measurement.  May be a linear function, i.e. y =
  %    A*x.
  % 
  %  KF_ic_state_func-function handle.  Takes information stored in
  %    KF_ic_state_func and turns it into initial conditions for the
  %    system.
  %
  %  KF_ic_paras_func-function handle.  Returns all parameters used in
  %    model (both the free ones that we are allowing to change and the set
  %    or static ones) as well as the corresponding standard deviations for
  %    the free parameters in a structure.  See 'ics_for_params_from_ppt'
  %    for details.
  %
  %  KF_noise_func-function handle, takes state vector and parameter value 
  %    ics and creates a noise structure with measurement, free parameter, 
  %    and state noise.  Also returns model and measurement error matrices.
  %
  %  KF_ic_state_addl_info-Additional information to feed to
  %    KF_ic_state_func function in parameter pairs, e.g.
  %    {'dx',2,'maxval',2}
  %
  %  KF_ic_paras_addl_info-cell of cells, e.g.:
  %    { {'kp',2,0.1,'free'},{'kpg',3},{'kq',[],[],'set'},...
  %      {parameter_name,parameter_value,...
  %       percent_parameter_stddev (optional),...
  %       parameter_free_or_set (optional)} }
  %    Allows us to adjust parameter info from default values.
  %    Note that if any of the slots are empty in the innermost cells, the
  %      default values for the particular parameter property are used.
  %      Also note, if the slot doesn't exist, again the default values are
  %      used.  Always the parameter name must appear in the first slot,
  %      the parameter value in the second slot, the parameter standard
  %      deviation in the third slot, and the choice of whether or not the
  %      parameter is to be free (and found via the KF) or set ('free' or
  %      'set') in the fourth slot.  Also, the standard deviation is a
  %      percentage/relative (not absolute) value.
  %  
  %  KF_noise_addl_info-Additional information in parameter pair form for
  %    noise adjustments from default values, e.g.
  %    {'model_noise',0.1,'kpg',0.1}.  Note the second arguments correspond
  %    to fractional standard deviations (relative to the mean parameter
  %    value) and that we can adjust the standard deviations corresponding
  %    to just one of our parameters...or all using 'kpg' vs 'param_noise'.
  %    Priority is always given to the individual vs the general so if we
  %    have {'param_noise',0.1,'kpg',0.2} or {'kpg',0.2,'param_noise',0.1},
  %    the standard deviation for 'kpg' will always be 20% mean value while
  %    the other std deviations will be 10% mean value.  Initial conditions
  %    can be adjusted by prepending 'ic', e.g. {'ic_state_noise',0.01}.  
  %    See KF_noise for details on possible changes.
  %  
  %  KF_time_info-again, parameter pairs of info correspond to time,
  %    e.g. {'texps',[1,2,3],'dtmax',4,'ts',[0:0.5:3]}
  %
  %  KF_options_cell-any additional options that we may want to change from
  %    their default values, in parameter pairs, e.g {'nm',0.1,...
  %    'em',0} (em is "experimental measurements" and the decision whether
  %    or not to use real or simulated data while nm is "noisy
  %    measurements" and the decision whether to add noise to the simulated
  %    data (it is never added to the real data))
  %
  %  problem_location-string giving the location of the files to be used
  %    during the Kalman filtering process
  
  %  Add the files to be used to the directory
  addpath(problem_location);

  %%  Assignment/reassignment of default values
  %  Initialize KF_options
  KF_options = get_options(KF_options_cell,KF_choice);
  
  %  Obtains a list of both set and unset parameters (params) and a list
  %  of just free parameters in state vector format
  [params,x0_params] = feval(KF_ic_params_func,KF_ic_params_addl_info{:});
  
  %  Obtain the initial conditions for the state of the system.  We have
  %  two possibilities
  %    1.  We want an initial state that corresponds to some experimental
  %      picture of the actual system.  In this case, inside the following
  %      program we not only extract from experiment the initial
  %      experimental data to be used in getting the original state of the
  %      system, we also, while we're at it, extract the rest of the
  %      experimental data as well.  This data is stored in exp_info.
  %    2.  We may want, for simplicity, some simple initial state
  %      corresponding to, e.g. a circle and not corresponding to
  %      experiment.  In this case we only find the initial state of the
  %      system and do not generate any of the exp_info.  In that case,
  %      exp_info = [] when handed back.
  %
  %  Note :  addl_info is a structure that carries around useful
  %    information.  Ideally it will contain info for the grid in 2-d/3-d
  %    simulations. It may be useful to allow it to carry other information
  %    too.
  [x0_mean,vs,exp_info,addl_info] = ...
    construct_state_vec_from_state_and_params(...
    KF_ic_state_func,x0_params,KF_ic_state_addl_info,KF_options);  
  
  %  This gets baseline levels for all the uncertainties (they may be
  %    adjusted along the way)
  noise_struc = KF_baseline_noise(params,KF_noise_addl_info{:});
  
  %  This retrieves simulated measurements for all time.  Note that our
  %  initial "guess" for the parameter values is stored in the state vector
  %  in x0_mean(vs.ns+1:vs.ns+vs.nf) while the exact values we desire to
  %  use for the simulated measurements are stored in params.free.val.
  %  That is why we make the call below the way that we do, to use the
  %  right parameter values for the simulated measurements.
  if isempty(exp_info)
    exp_info = feval(KF_exp_func, ...
      [x0_mean(1:vs.ns);cell2mat({params.free.val})'], ...
      params, KF_options, noise_struc, addl_info, exp_info, ...
      KF_meas_func, KF_meas_addl_info, vs, KF_time_info{:});
  end
  
  %  This gets the initial mean and additional useful info (as described
  %    above)
  
  %  These are the times at we we have experimental measurements...also the
  %    times at which we want assimilate the data/use the kalman filter
  %    step
  ts = exp_info.expts;
  %  Note we assume the amount of data available does not change over time
  %    (may not always be true--it's easy to change the code when it is not
  %    the case).
  vs.m = size(exp_info.meas_vecs,1);
  
  %  This gets initial covariances, baseline levels
  initial_baseline_noise = [noise_struc.stddev.ini_state*ones(1,vs.ns),...
    noise_struc.stddev.ini_params];
  Pxx_0 = diag(initial_baseline_noise.^2);
  
  %  This gets model covariances, baseline levels
  comp_baseline_noise = [noise_struc.stddev.model_state*ones(1,vs.ns),...
    noise_struc.stddev.model_params];
  Pxx_model = diag(comp_baseline_noise.^2);
  
  %  This gets measurement covariances, baseline levels
  meas_baseline_noise = noise_struc.stddev.meas*ones(1,vs.m);
  Pyy_meas = diag(meas_baseline_noise.^2);
  
  %  Generate random or organized ensembles as need be
  switch KF_options.KF_type
    case 'En_KF'
      xa_en = bsxfun(@plus,x0_mean,diag(initial_baseline_noise)*...
        randn(vs.n,KF_options.q));
    case 'SC_KF'
      [weights,organized_normalized_noise] = Jared_Sparse(vs.n);
      xa_en = bsxfun(@plus,x0_mean,diag(initial_baseline_noise)*...
        organized_normalized_noise);
    case 'SC_SVD_KF'
      [weights,organized_normalized_noise] = Jared_Sparse(vs.n);
      Us = sparse([],[],[],vs.n,vs.n,0);
      Vs = Us;
      Ss = Vs;
      sigsn = min(KF_options.M2,vs.n);
      [U,S,V] = svds(Pxx_model,sigsn);
      Us(:,1:sigsn) = U; Vs(:,1:sigsn) = V; Ss(1:sigsn,1:sigsn) = S;
      Pxx_model_abbrev = Us*Ss*Vs';
      Pxx_model_abbrev_sqrt = Us*sqrt(Ss)*Vs';
      xa_en = bsxfun(@plus,x0_mean,Pxx_model_abbrev_sqrt*...
        organized_normalized_noise);
    case 'En_SVD_KF'
      Us = sparse([],[],[],vs.n,vs.n,0);
      Vs = Us;
      Ss = Vs;
      sigsn = min(KF_options.M2,vs.n);
      [U,S,V] = svds(Pxx_model,sigsn);
      Us(:,1:sigsn) = U; Vs(:,1:sigsn) = V; Ss(1:sigsn,1:sigsn) = S;
      Pxx_model_abbrev = Us*Ss*Vs';
      Pxx_model_abbrev_sqrt = Us*sqrt(Ss)*Vs';
      xa_en = bsxfun(@plus,x0_mean,Pxx_model_abbrev_sqrt*...
        randn(vs.n,KF_options.q));
  end
  
  %  Initialize important parts of the KF algorithm
  Pxxa = Pxx_0;
  xa_mean = x0_mean;
  
  %  Store somewhere
  Pxxas = zeros([size(Pxxa),length(ts)]);
  Pxxas(:,:,1) = Pxxa;
  xa_means = zeros(length(xa_mean),length(ts));
  xa_means(:,1) = xa_mean;
  
  for tc = 2:length(ts)
    %  Take one forward euler step.  The function should also return (for
    %  linear kalman filter/extended kalman filter) the jacobian
    %  corresponding to the model.
    dt = ts(tc)-ts(tc-1);
    
    %  Forecast step including getting the corresponding measurements and,
    %  if we are using the linear or extended KF, the jacobian of the
    %  measurement function
    %  For ensemble KF also calculate the means and error functions
    switch KF_options.KF_type
      case 'KF'
        [xf_mean,A] = feval(KF_model_func,xa_mean,dt,params,addl_info,vs);
        [yf_mean,H] = feval(KF_meas_func,xf_mean,vs,KF_meas_addl_info);
      case 'En_KF'
        [xf_en,trash] = feval(KF_model_func,xa_en,dt,params,addl_info,vs);
        xf_en =  bsxfun(@plus, xf_en, ...
          bsxfun(@times,comp_baseline_noise',randn(vs.n,KF_options.q)));
        xf_mean = mean(xf_en,2);
        [yf_en,trash] = feval(KF_meas_func,xf_en,vs,KF_meas_addl_info);
        yf_mean = mean(yf_en,2);
      case 'En_SVD_KF'
        [xf_en,trash] = feval(KF_model_func,xa_en,dt,params,addl_info,vs);
        xf_mean = mean(xf_en,2);
        [yf_en,trash] = feval(KF_meas_func,xf_en,vs,KF_meas_addl_info);
        yf_mean = mean(yf_en,2);
      case {'SC_KF','SC_SVD_KF'}
        [xf_en,A] = feval(KF_model_func,xa_en,dt,params,addl_info,vs);
        xf_mean = sum(xf_en*diag(weights),2);
        [yf_en,H] = feval(KF_meas_func,xf_en,vs,KF_meas_addl_info);
        yf_mean = sum(yf_en*diag(weights),2);
    end
    
    %  Calculate new covariances including the following steps
    %  1)  Calculate the new covariance due to using the model over that one time
    %  step
    %  2)  Calculate the covariance between the state vector/model and the
    %  measurements according to the model
    %  3)  Calculate the covariance between the state vector/model and the
    %  measurements according to the model
    switch KF_options.KF_type
      case 'KF'
        Pxxf = A*Pxxa*A'+Pxx_model;
        Pxyf = Pxxf*H';
        Pyyf = H*Pxxf*H'+Pyy_meas;
      case 'En_KF'
        Exf = bsxfun(@minus,xf_en,xf_mean);
        Eyf = bsxfun(@minus,yf_en,yf_mean);
        Pxxf = (1/(KF_options.q-1))*Exf*transpose(Exf);
        Pxyf = (1/(KF_options.q-1))*Exf*transpose(Eyf);
        Pyyf = (1/(KF_options.q-1))*Eyf*transpose(Eyf);
      case 'En_SVD_KF'
        Exf = bsxfun(@minus,xf_en,xf_mean);
        Eyf = bsxfun(@minus,yf_en,yf_mean);
        Pxxf = (1/(KF_options.q-1))*Exf*transpose(Exf)+Pxx_model;
        Pxyf = (1/(KF_options.q-1))*Exf*transpose(Eyf);
        Pyyf = (1/(KF_options.q-1))*Eyf*transpose(Eyf)+Pyy_meas;
      case {'SC_KF','SC_SVD_KF'}
        Exf = bsxfun(@minus,xf_en,xf_mean);
        Eyf = bsxfun(@minus,yf_en,yf_mean);
        Pxxf = Exf*diag(weights)*(Exf')+Pxx_model;
        Pxyf = Exf*diag(weights)*(Eyf');
        Pyyf = Eyf*diag(weights)*(Eyf')+Pyy_meas;
%         Pyyf = Pyyf';
    end
 
    %  Calculate the Kalman gain...the magic blending factor
    Kn = Pxyf/(Pyyf);
    
    %  Calculate the new adjusted state vector...the new mean that is the
    %  optimal blend of the old mean according to the model and the
    %  experimental measurements.  For the ensemble KF, do the correction
    %  for each member of the ensemble
    switch KF_options.KF_type
      case {'KF','SC_KF','En_SVD_KF','SC_SVD_KF'}
        xa_mean = xf_mean+Kn*(exp_info.meas_vecs(:,tc)-yf_mean);
      case {'En_KF'}
        xa_en = xf_en+...
          Kn*(...
          bsxfun(@plus,exp_info.meas_vecs(:,tc),...
          diag(meas_baseline_noise)*randn(vs.m,KF_options.q))-yf_en);
        xa_mean = mean(xa_en,2);
    end
    
    %  Different formulas for the new covariance associated with the new
    %  adjusted state vector estimate (various papers use different formulas
    %   Pxxa = (eye(size(Pxxf))-Kn*H)*Pxxf;
    %   Pxxa = Pxxf-Kn*Pxdf';%(eye(size(Pxxf))-Kn*H)*Pxxf;
    %  Gillijns
    switch KF_options.KF_type
      case 'KF'
        Pxxa = Pxxf-Pxyf*Kn';%Kn*H*Pxxf;
      case 'En_KF'
        Exa = bsxfun(@minus,xa_en,xa_mean);
        Pxxa = (1/(KF_options.q-1))*Exa*transpose(Exa);
      case 'En_SVD_KF'
        Pxxa = Pxxf-Pxyf*transpose(Kn);
        Pxxa = (Pxxa+Pxxa')./2;
        [U,S,V] = svds(Pxxa,sigsn);
        Us(:,1:sigsn) = U; Vs(:,1:sigsn) = V; Ss(1:sigsn,1:sigsn) = S;
        Pxxa_abbrev_sqrt = Us*sqrt(Ss)*Vs';
        xa_en = bsxfun(@plus,xa_mean,Pxxa_abbrev_sqrt*...
          randn(vs.n,KF_options.q));
      case 'SC_KF'
        Pxxa = Pxxf-Pxyf*Kn';
%         Pxxa = (Pxxa+Pxxa')./2;
        xa_en = bsxfun(@plus,xa_mean,...
          diag(sqrt(diag(Pxxa)))*organized_normalized_noise);
%         fprintf('%g ',Pxxa-Pxxa');
%         fprintf('\n');
      case 'SC_SVD_KF'
        Pxxa = Pxxf-Pxyf*transpose(Kn);
        Pxxa = (Pxxa+Pxxa')./2;
        [U,S,V] = svds(Pxxa,sigsn);
        Us(:,1:sigsn) = U; Vs(:,1:sigsn) = V; Ss(1:sigsn,1:sigsn) = S;
        Pxxa_abbrev_sqrt = Us*sqrt(Ss)*Vs';
        xa_en = bsxfun(@plus,xa_mean,...
          Pxxa_abbrev_sqrt*organized_normalized_noise);
%         fprintf('%g ',diag(full(Ss)));
%         fprintf('\n');
        [U,S,V] = svd(Pxxa);
%         fprintf('%g ',diag(full(S)));
%         fprintf('\n');
    end
      
    %  Store all the vectors and matrices found for later use
    xf_means(:,tc) = xf_mean;
    xa_means(:,tc) = xa_mean;
    Pxxfs(:,:,tc) = Pxxf;
    Pxxas(:,:,tc) = Pxxa;
    
  end
  
  %  Plots
  if KF_options.plot_all
    
    close(figure(KF_options.fig_nums(1))); figure(KF_options.fig_nums(1));
    %  Plot all states
    subplot(3,max(2,vs.nf),1)
    plot(ts,xa_means(1:vs.ns,:),'b',...
      ts,exp_info.meas_vecs,'r')
    title('Blue-Computational.  Red-Simulated/Experimental');
    
    %  Plot differences in states (should get closer to zero as time goes
    %  on).
    if size(exp_info.meas_vecs,1) == vs.ns
      subplot(3,max(2,vs.nf),2)
      plot(ts,xa_means(1:vs.ns,:)-exp_info.meas_vecs(1:vs.ns,:),'b');
      title('Model measurements-Experimental Measurements');
    end
    
    %  Plot parameters
    for fpc = 1:vs.nf
      subplot(3,max(2,vs.nf),max(2,vs.nf)+fpc)
      plot(ts,xa_means(vs.ns+fpc,:),'b');
      title(params.free(fpc).name);
    end
    
    %  Plot parameters vs exact (only good for simulated data)
    if ~KF_options.em
      for fpc = 1:vs.nf
        subplot(3,max(2,vs.nf),2*max(2,vs.nf)+fpc)
        plot(ts,abs(xa_means(vs.ns+fpc,:)-params.free(fpc).val)./...
          params.free(fpc).val,'b');
        set(gca,'YScale','log');
        title(['Relative error in KF estimate for ',...
          params.free(fpc).name]);
      end
    end
  end
  
  %  Remove the files from the path after been done
  rmpath(problem_location);
  
end

function [x0,vs,exp_info,addl_info] = ...
  construct_state_vec_from_state_and_params(...
  KF_ic_state_func,x0_params,KF_ic_state_addl_info,KF_options)
  
  exp_info = [];
  %  Obtains the initial state vector for the system as well as other
  %  additional info.  The addl_info structure can contain anything but
  %  most importantly it usually contains grid information.  As noted
  %  earlier, if we want to use experimental data, this program will
  %  automatically use the experimental data to obtain the initial
  %  condition and return all experimental data in exp_info.  Otherwise,
  %  exp_info will return empty to be filled later with simulated data.
  [x0_state,addl_info,exp_info] = feval(KF_ic_state_func,exp_info,...
    KF_options,KF_ic_state_addl_info{:});
  
  %  Make the state vector complete with state values and parameter values
  x0 = [x0_state;x0_params];
  
  %  Store information about the vector size for the state.
  %  number of free parameters
  vs.nf = length(x0_params);
  
  %  number of state values
  vs.ns = length(x0_state);
  
  %  number of total free values
  vs.n = vs.nf+vs.ns;
end


function KF_options = get_options(KF_options_cell,KF_choice)

  %  Use experimental measurements option
  KF_options.em = false;
  
  %  Use noisy simulated measurements
  KF_options.nm = false;
  
  %  Which choice we want to use for KF, stochastic collocation,
  %  Karhunen-Loeve, etc.
  KF_options.KF_type = KF_choice;
  
  %  Whether or not we adjust the standard deviations for the state vector
  %  in order to guarantee all ensemble members have only positive values
  KF_options.positivity = false;
  
  %  Whether or not to plot state vector values and uncertainties at the
  %  end
  KF_options.plot_all = true;
  
  %  Figure number to use when plotting
  KF_options.fig_nums = [1];
  
  %  Number of ensemble members (for ensemble kalman filters)
  KF_options.q = 1000;
  
  %  Number of degrees of freedom (for KL or svd kalman filters), i.e. the
  %  number of eigendirections that we decide to include in our noise
  KF_options.M2 = 25;
%   
%   %  Whether or not to bump the initial condition for the model
%   KF_options.
  
  %  We assume the user hands in temp_KF_options in a cell...this loop
  %  replaces the above default values with user-defined values
  for tKoc = 1:2:length(KF_options_cell)
    KF_options.(KF_options_cell{tKoc}) = KF_options_cell{tKoc+1};
  end
  
  %  Set adding noise to "simulated measurements" to zero if we are using
  %  real data as real data already has noise in it.
  KF_options.nm = KF_options.nm & (~KF_options.em);
  
end

function noise_struc = KF_baseline_noise(params, varargin)

  %  Important note, uncertainties associated with parameters are dealt
  %    with when parameters are assigned in ics_for_params_from_ppt.
  
  %% Default values for noise
  %  noise_frac is nice to have around as you can dial down or up all
  %    noises at the exact same rate/time
  noise_frac = 1;
  model_state_noise = 0.003;
  ini_state_noise = 0.003;
  meas_noise = 0.003;
  
  %%  Readjustment of those default values
  for vac = 1:2:length(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  
  %%  Extract uncertainties associated with parameters
  %  Note: the assumptions listed below aren't always true.
  param_noise = [];
  for fpc = 1:length(params.free)
    param_noise(fpc) = params.(params.free(fpc).name).stddev;
  end
  %  Assume the initial parameter uncertainty is given by the parameter
  %    noises above
  ini_param_noise = param_noise;
  
  %  Assume the uncertainty associated with using the model to adjust the
  %    parameter values is also given by the above parameter noises
  model_param_noise = param_noise;
  
  %%  Store this information in baseline covariance matrix
  %  Note that these are "baseline" standard deviations/uncertainties.
  %    Their values may change during the actual Kalman filtering process.

  %  Std deviations for computational and data noise
  noise_struc.stddev.model_state = noise_frac*model_state_noise;
  noise_struc.stddev.ini_state = noise_frac*ini_state_noise;
  noise_struc.stddev.meas = noise_frac*meas_noise;
  noise_struc.stddev.model_params = noise_frac*model_param_noise;
  noise_struc.stddev.ini_params = noise_frac*ini_param_noise;

  %  When these quantities are nonzero, the noise will decay to zero as
  %  time goes on.  In particular, noise is proportional to
  %  exp(-decrate*t).  This makes sense to do for parameters...as time goes
  %  on we should get closer and closer to the correct parameter estimates.
  %  As such our uncertainties associated with the parameter values
  %  decrease with time.  While it doesn't generally make sense to decrease
  %  the uncertainties associated with the model and/or measurements as
  %  time goes on, there may be certain situations where decreasing those
  %  uncertainties does make sense.
  noise_struc.decrate.model_state = 0;
  noise_struc.decrate.meas = 0;
  noise_struc.decrate.model_params = 0;

end