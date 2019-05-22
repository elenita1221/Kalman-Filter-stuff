function [grid_weight, grid_point] = Jared_Sparse (dim_num)

  compare_on = 0;
  
  w = (1/6)*ones(1,2*dim_num+1);
  w(dim_num+1) = 1-dim_num/3;
  n = [-sqrt(3)*speye(dim_num);sparse(zeros(1,dim_num));...
    sqrt(3)*flipud(speye(dim_num))];
  
  if compare_on == 1
    ad = cd; cd('C:\NEC\Sparse_Grid_2'); 
    [n2 w2] = nwspgr('KPN',dim_num,2); cd(ad);
    fprintf('Comparison: %g\n',max(max(abs(w(:)-w2(:))),max(n(:)-n2(:))));
  end
  
  grid_weight = w;
  grid_point = n';

end