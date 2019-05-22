function ecs = ec_model(params, free_param_vec,free_param_names,...
  ecini,ets,g,varargin)

%  ecs = ec_model(param_vec,ecini,ts,g)
%
%  free_param_vec-vector of free parameters, typically param_vec = [D,kp], 
%    just two though it can be adjusted as need be so long as you adjust
%    the program accordingly, in particular the subprogram advance_ecs
%    would need to be changed.
%  free_param_names-names of the corresponding free parameters listed in
%    same order as their values appear in the free_param_vec
%  ecini-initial epithelial cell density in the grid
%  ets-experimental times at which we want to know the epithelial cell 
%    density in the grid
%  g-grid info including
%    g.nx-number of cells in x-direction
%    g.ny-number of cells in y-direction
%    g.dx-size of cells in x-direction (we assume uniform spacing)
%    g.dy-size of cells in y-direction (again uniform spacing assumed)
%  Note:  The info inside g is later augmented with the appropriate time
%    step to use during each step in the numerical process (changes with
%    time).

%%  Set default parameter values
%  D-diffusion coefficient, corresponds to diffusion rate of cells
%  kp-growth or proliferation rate of epithelial cells
%  hc-hill coefficient, power in the hill term in the diffusion
%  Note:  To make the hill coefficient a free parameter, you need to swap
%  the commented and uncommented lines below for hc.  You would also need
%  to make the right changes in the MCMC program that called this program.
%  Note:  Currently the parameter values are handed in
%   params.D.val = 3e-6;
%   params.kp.val = 2;
%   params.hc.val = 2;
  
  %%  If the above parameters are free, replace their values
  %  by those that were handed into the program from the KF or from MCMC
  [xM,yM] = meshgrid(g.xcs,g.ycs);
  for fpc = 1:length(free_param_names)
    params.(free_param_names{fpc}).val = free_param_vec(fpc);
  end
  
 % params.D.val= ones(g.ny,g.nx)+xM.^2+yM.^2;
%for i=1:5 
%     params.D.val=params.D.val+cos(i*pi*xM/0.1).*cos(i*pi*yM/0.07).*(1/2);
%end
%params.D.val=params.D.val*(1e-3);

  %%  Consts for good numerical solving and other options
  ecm.safety_net =0.5;
  ecm.plot_soln = 1;
  
  for vac = 1:2:length(varargin)
    ecm.(varargin{vac}) = varargin{vac+1};
  end

  %%  Get maximum dt allowable for stability
  %  This includes a safety buffer fraction of 0.9.
  temp_quant = g.dx^2*g.dy^2/(params.D.val*(g.dx^2+g.dy^2));
  maxdt = ecm.safety_net*temp_quant;
  
  ecs = zeros([size(ecini),length(ets)]);
  ecs(:,:,1) = ecini;
  my_plot(ecm.plot_soln,ecini)

  for etc = 1:length(ets)-1

    exp_dt = ets(etc+1)-ets(etc);
    nsteps = ceil(exp_dt/maxdt);
    comp_dt = exp_dt/nsteps;
    g.dt = comp_dt;

    ecs(:,:,etc+1) = ecs(:,:,etc);
    for ctc = 1:nsteps
      ecs(:,:,etc+1) = advance_ecs(ecs(:,:,etc+1),g,params);
    end

    my_plot(ecm.plot_soln,ecs(:,:,etc+1))

  end
  
end

function my_plot(plot_soln,ecs)
  if plot_soln
    surf(ecs);
    zlim([0,1]);
    hold on;
    contour(ecs,[0.5 0.5]);
    hold off;
    pause;
  end
end

function ecs = advance_ecs(ecs,g,params)
  
 grad_x = [(ecs(1,:)-0)./(g.dx/2);...
    diff(ecs,1,1)./g.dx;...
    (0-ecs(end,:))./(g.dx/2)];

 grad_y = [(ecs(:,1)-0)./(g.dy/2),...
    diff(ecs,1,2)./g.dy,...
    (0-ecs(:,end))./(g.dy/2)];
  
ecs = ecs+g.dt*(params.D.val*(diff(grad_x,1,1)./g.dx+...
    diff(grad_y,1,2)./g.dy)+0*ones(g.ny,g.nx));

end
