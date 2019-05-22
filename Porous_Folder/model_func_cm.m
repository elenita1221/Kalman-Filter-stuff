function [xnext,A] = model_func_cm(xprev,dt,params,addl_info,vs)
  %  x-state vector
  %  dt-time step for the model
  %  params-structure containing parameter information for set parameters.
  %    In this function no parameters are set so the empty set is handed in
  %    for params and it is never used in this program.
  %  addl_info-structure containing additional useful information including
  %    grid information
  %  Note:  Free parameters are handed in at the end of xprev
  A = [];
  xnext = xprev;
  nx = addl_info.g.nx;
  ny = addl_info.g.ny;
  g = addl_info.g;
  for enc = 1:length(xprev)
    temp = xprev(:,enc);
    temp_ecs = ec_model(params,temp(vs.ns+1:vs.n),...
      {params.free.name},reshape(temp(1:vs.ns),nx,ny),...
      [0,dt],g,'plot_soln',0);
    temp_ecs = temp_ecs(:,:,end);
    xnext(:,enc) = [temp_ecs(:);temp(vs.ns+1:vs.n)];
  end

end