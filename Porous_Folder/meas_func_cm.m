function [meas,H] = meas_func_cm(x,vs,meas_addl_info)

  %  Note x can be a vector or a matrix of column vectors where each column
  %  is the corresponding state that we want to find the measurement of.
  meas = x(1:vs.ns,:);
  H = eye(vs.ns,vs.n);
end