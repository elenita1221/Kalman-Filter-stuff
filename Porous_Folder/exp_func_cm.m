function exp_info = exp_func_cm(x0_mean, params, KF_options, noise_struc,...
  addl_info, exp_info, meas_func, meas_addl_info, vs, varargin)
  %  varargin-list of options specific to this program.  The default values
  %    for various parameters are given below.  All of these default values
  %    may be changed by handing in parameter pairs in varargin.  e.g.
  %    simulated_f_from_ppt('expts',[0,1,2,3]) will change the default values
  %    for "ts" below to be [0,1,2,3].

  %%  Default values
  %  Experimental times.  Note these times correspond only to the simulated
  %  data as these default values are replaced later by the actual data
  %  from experiment if we are not using simulated data.
  exp_info.expts = linspace(0,3.75,16);
  
  %  Sometimes we bump the simulated measurements by some noise.  This is
  %  the standard deviation of the corresponding noise bump and if it is 0,
  %  then we use "perfect simulated measurements" that are not bumped at
  %  all.
  exp_info.bump_size = 0.1;
  
  %  A label that this is simulated, not real, measurements
  exp_info.rm = 0;
  
  %  First exact state vector (not a measurement, includes params)
  exp_info.xe0 = x0_mean';
  
  %  Function that changes state vector into a measurement vector
  exp_info.meas_func = meas_func;
  
  %  Function that advances solution from one time step to the next, given
  %  the state vector, time step, and parameter structure
  exp_info.model_func = @f_from_ppt;
  
  %%  Reset default values according to desired values
  for vac = 1:2:length(varargin)
    exp_info.(varargin{vac}) = varargin{vac+1};
  end
  
  %%  Reset default values to agree with experiment
  if KF_options.em
    exp_info.expts = exp_info.TimeExperiments;
  end
  
  %%  Run an "exact" simulation in order to produce measurements
 
  %  Initialize a matrix whose columns correspond to the simulated
  %  state vector values
  state_vec_matrix = zeros(length(exp_info.xe0),length(exp_info.expts));
  
  %  Initialize the first column of the matrix, corresponds to the initial
  %  condition
  state_vec_matrix(:,1) = exp_info.xe0;

  %  If experimental measurements = false, then go ahead and make simulated
  %  measurements.
  if ~KF_options.em
    temp_ecs = ec_model(params,exp_info.xe0(vs.ns+1:vs.n),...
      {params.free.name},reshape(exp_info.xe0(1:vs.ns),addl_info.g.nx,...
      addl_info.g.ny),exp_info.expts,addl_info.g,'plot_soln',0);
    temp_ecs = permute(temp_ecs,[3,1,2]);
    
  %  Otherwise, fetch the experimental measurements (loaded in earlier
  %  somewhere else)
  else
    temp_ecs = exp_info.ExpData;
  end
  
  %  In both cases, we reshape this experimental/simulated measurements
  %  into state(/measurement) vector format
  state_vec_matrix(1:vs.ns,:) = reshape(temp_ecs,length(exp_info.expts),...
    vs.ns)';
  
  %  
  for fpc = 1:vs.nf
    state_vec_matrix(vs.ns+fpc,:) = exp_info.xe0(vs.ns+fpc);
  end
  
  exp_info.exact_state = state_vec_matrix;
  
  %  Initialize a matrix whose columns correspond to the simulated
  %  measurements
  exp_info.meas_vecs = feval(exp_info.meas_func,state_vec_matrix,vs,...
    meas_addl_info);
  
  %  Bump measurements to simulate imperfect measurements in reality
  exp_info.meas_vecs = exp_info.meas_vecs+exp_info.bump_size*...
    noise_struc.stddev.meas*KF_options.nm*...
    normrnd(zeros(size(exp_info.meas_vecs)),ones(size(exp_info.meas_vecs)));
  
end