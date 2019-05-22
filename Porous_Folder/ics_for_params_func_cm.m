function [params,x0_params] = ics_for_params_func_cm(varargin)

  %%  Some default values/options for the program
  %   The measurement noise on the parameters is assumed to be gaussian
  %   distribution with a standard deviation = stddev_frac*original value
  %   of the parameter in question (see calculations below).  When
  %   stddev_frac = 0.1, the standard deviation is equal to 10% of the
  %   original value of the parameter
  stddev_frac =sqrt(0.1);
  
  %%  Readjust those default values
  varargin_new = {};
  for vac = 1:2:length(varargin)
    switch varargin{vac}
      case 'stddev_frac'
        stddev_frac = varargin{vac+1};
      otherwise
        varargin_new = {varargin_new{:},varargin{vac:vac+1}};
    end
  end

  %%  Define default values for parameters
  %  type-> whether or not a parameter is free.  Parameters can also be
  %    "set" and not allowed to roam.
  %  val-> The exact value of the parameter used when producing simulated
  %    measurements (never used if real measurements are used)
  %  initial_guess-> Guess for the exact value of the measurement (making
  %    right guess can be crucial when using real measurements)
  %  stddev-> The estimated uncertainty associated with the parameter
 %[xM,yM] = meshgrid(g.xcs,g.ycs);
 % params.D.val=0.5*ones(g.ny,g.nx);
 % params.D=params.D.val.*(1e-3);
  params.D.type = 'free';
  params.D.val = 3e-3;
  params.D.initial_guess = 1.5e-3;
 % params.D.stddev = stddev_frac*params.D.initial_guess;
  params.D.stddev = stddev_frac*params.D.val;
  
 % params.kp.type = 'free';
 % params.kp.val = 1;
 % params.kp.initial_guess = 0.5;
 % params.kp.stddev = stddev_frac*params.kp.val;
  
  params.hc.type = 'set';
  params.hc.val = 2;
%   params.hc.initial_guess = 2;  Not needed when parameter is set.
%   params.hc.stddev = stddev_frac*params.hc.val;
  
  %%  Readjust default values
  %  Note varargin is in parameter pair form again,
  %  ics_for_params_from_ppt('a.val',2) would readjust the variable a used
  %  in exact simulated measurements to be 2 while
  %  ics_for_params_from_ppt('a.initial_guess',1) would change the initial
  %  guess for the variable a to be 1.
  for vac = 1:2:length(varargin_new)
    tmp_ind = findstr('.',varargin_new{vac});
    if isempty(tmp_ind)
      error(['When readjusting a parameter, you must specify if you ',...
        'wish to readjust the value used during simulated ',...
        'measurements (.val) or the value used as the initial guess for',...
        'the parameter value (.initial_guess)']);
    else
      param_name = varargin_new{vac}(1:tmp_ind-1);
      param_val_or_guess = varargin_new{vac}(tmp_ind+1:end);
      params.(param_name).(param_val_or_guess) = varargin_new{vac+1};
      if strcmpi(param_val_or_guess,'val')
        params.(param_name).stddev = stddev_frac*varargin_new{vac+1};
      end
    end
  end
  
  %%  Store default values in an easier to use location  
  %  So that it's easy to stick these things into and take them out of
  %  state vectors we arrange them inside of their own structures
  fn = fieldnames(params);
  
  %  Counter for free parameters
  fc = 0;
  
  %  Loop through and store the free parameters and corresponding info in
  %  the "free" field
  for fnc = 1:length(fn)
    if strcmpi('free',params.(fn{fnc}).type) == 1
      fc = fc+1;
      params.free(fc).name = fn{fnc};
      params.free(fc).val = params.(fn{fnc}).val;
      params.free(fc).stddev = params.(fn{fnc}).stddev;
    end
  end
  if fc == 0, params.free = []; end
  
  %  Do a similar thing for set parameters.  Note, since the parameters are
  %  set, they have no associated stddev.

  %  Counter for set parameters
  sc = 0;
  
  %  Loop through and store the set parameters and corresponding info in
  %  the "set" field
  for fnc = 1:length(fn)
    if strcmpi('set',params.(fn{fnc}).type) == 1
      sc = sc+1;
      params.set(sc).name = fn{fnc};
      params.set(sc).val = params.(fn{fnc}).val;
    end
  end
  if sc == 0, params.set = []; end
  
  %  Store all initial guesses for the free parameters in state vector
  %  format
  for fpc = 1:length(params.free)
    x0_params(fpc) = params.(params.free(fpc).name).initial_guess;
  end
  x0_params = reshape(x0_params,length(x0_params),1);
  
end