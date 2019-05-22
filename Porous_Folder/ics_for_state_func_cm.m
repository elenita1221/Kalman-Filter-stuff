function [x0_state,addl_info,exp_info] = ics_for_state_func_cm(exp_info,...
  KF_options, varargin)

  %  In this particular program, varargin will never be used, but I just
  %  wanted to include them as they may be included later and may include
  %  things like directories where experimental data is found and/or
  %  programs that take experimental data and convert them somehow into an
  %  initial state vector.
  %
  %  params gives the initial guesses for the parameter values (decided
  %  earlier in the KF_ic_paras_func function)
  %
  %  Note this is not the exact vector that corresponds to the simulated
  %  measurements.  Rather, it is off by a bit, in the parameter estimate
  %  area in particular.  One could also perturb the other state vector
  %  entries to see what happens.
  addl_info = [];
  
  %%  Some default values that can easily be readjusted
  %   Note all the other grid information pretty much depends on these
  %   guys.  Grid info:
  addl_info.g.nx = 10;
  addl_info.g.ny = 10;
  addl_info.g.xext = 0.1;
  addl_info.g.yext_div_xext = 0.07/0.1;
  
  %  We may choose different hole types to investigate as well as different
  %  initial depths and areas 
  %  type->options = 'circle' and 'irregular' where the latter correspodns
  %  to experiment 
  %  fracarea->initial fractional area of the injury.  Note for circular
  %  hole and rectangular domains, the area < approximately 60% to make
  %  sure the circular injury does not hit the domain boundaries.  Also we
  %  do not yet readjust the area of the original wound for irregular holes
  %  severity->Usually 1 in the injured region
  addl_info.ht.type = 'circle'
  addl_info.ht.fracarea = 0.2;
  addl_info.ht.severity = 1;
  
  %%  Readjust default values
  for vac = 1:2:length(varargin)
    tmp_ind = findstr('.',varargin{vac});
    if isempty(tmp_ind)
      error(['When readjusting a addl_info field, you must specify',...
        'the substructure you wish to adjust.  E.g. to readjust',...
        'addl_info.g.nx, you must hand in ''g.nx'', not just ''nx''.']);
    else
      substruc_name = varargin{vac}(1:tmp_ind-1);
      field_name = varargin{vac}(tmp_ind+1:end);
      if isfield(addl_info.(substruc_name),field_name)
        addl_info.(substruc_name).(field_name) = varargin{vac+1};
      else
        error('You entered an invalid field for the addl_info structure!');
      end
    end
  end
  
  %%  Readjust default values based on experiment
  %  If experimental measurements are being used we force the nx, ny, xext,
  %  and yext_div_xext to match up with what is recorded in the data file.
  %  We also print a warning to that extent.  For the time being, we assume
  %  the default values set above are actually more accurate than the ones
  %  recorded (currently true I believe).  For that reason everything below
  %  is commented out for now.
  if KF_options.em == 1
    exp_info = load('KalmanDataNEC_nx10_ny10.mat');
    exp_info.expts = exp_info.TimeExperiments;
    for tc = 1:length(exp_info.expts)
      ec_temp = exp_info.ExpData(tc,:,:);
      exp_info.meas_vecs(:,tc) = ec_temp(:);
    end
    %  Only select some of the experimental times (gives more monotonic
    %  results)
    tinds = 1:3:48;
    exp_info.expts = exp_info.expts(tinds);
    exp_info.meas_vecs = exp_info.meas_vecs(:,tinds);
    
%     if exp_info.XCells ~= addl_info.g.nx
%       addl_info.g.nx = exp_info.XCells;
%       warning(['Requested number of cells in the x-direction and ',...
%         'number of cells in the x-direction in the experimental files',...
%         'do not agree!  Using experimental information.']);
%     end
%     if exp_info.YCells ~= addl_info.g.ny
%       addl_info.g.ny = exp_info.YCells;
%       warning(['Requested number of cells in the y-direction and ',...
%         'number of cells in the y-direction in the experimental files',...
%         'do not agree!  Using experimental information.']);
%     end
%     if abs(exp_info.exp_ext(1)-addl_info.g.xext)>1e-10
%       addl_info.g.xext = exp_info.exp_ext(1);
%       warning(['Requested domain extent in the x-direction and ',...
%         'domain extent in the x-direction in the experimental files',...
%         'do not agree!  Using experimental information.']);
%     end
%     if abs(exp_info.exp_ext(2)-...
%         addl_info.g.xext*addl_info.g.yext_div_xext)>1e-10
%       addl_info.g.yext_div_xext = exp_info.exp_ext(2)/exp_info.exp_ext(1);
%       warning(['Requested domain extent in the y-direction and ',...
%         'domain extent in the x-direction in the experimental files',...
%         'do not agree!  Using experimental information.']);
%     end
  end
  
  %%  Now set up other grid information used later in calculations
  addl_info.g.yext = addl_info.g.yext_div_xext*addl_info.g.xext;
  addl_info.g.xl = -addl_info.g.xext/2;
  addl_info.g.xr = addl_info.g.xext/2;
  addl_info.g.dx = (addl_info.g.xr-addl_info.g.xl)./addl_info.g.nx;
  addl_info.g.yl = -addl_info.g.yext/2;
  addl_info.g.yr = addl_info.g.yext/2;
  addl_info.g.dy = (addl_info.g.yr-addl_info.g.yl)./addl_info.g.ny;
  
  addl_info.g.xes = ...
    linspace(addl_info.g.xl,addl_info.g.xr,addl_info.g.nx+1);
  addl_info.g.yes = ...
    linspace(addl_info.g.yl,addl_info.g.yr,addl_info.g.ny+1);
  addl_info.g.xcs = (addl_info.g.xes(1:end-1)+addl_info.g.xes(2:end))./2;
  addl_info.g.ycs = (addl_info.g.yes(1:end-1)+addl_info.g.yes(2:end))./2;
  
  
  
  %%  Use this to define ics
  %  Note we assume that if we use experimental data that we obtain our
  %  initial conditions by using the experimental data.
  if KF_options.em == 1
    ecini = squeeze(exp_info.ExpData(1,:,:));
  else
    switch addl_info.ht.type
      case 'circle'
        [xM,yM] = meshgrid(addl_info.g.xcs,addl_info.g.ycs);
         ecini = double((xM.^2+yM.^2)<0.03^2);
        % ecini = double((xM.^2+yM.^2)>...
        %  (addl_info.ht.fracarea*addl_info.g.xext*addl_info.g.yext/pi));
      case 'irregular'
        temp_data = load('KalmanDataNEC_nx10_ny10.mat');
        ecini = temp_data.ExpData(1,:,:);
    end
  end
  x0_state = ecini(:);

end