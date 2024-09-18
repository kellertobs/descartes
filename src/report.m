% get residual of phase field update
resnorm = norm(upd_zp(:)) / (norm(zp(:)) + eps);

% Reset reference residual if necessary
if iter == 1 || resnorm > resnorm0
    resnorm0 = resnorm + eps;
end

% Check for solver failure (NaN in residuals)
if isnan(resnorm)
    error('!!! Solver failed with NaN: end run !!!');
end

% Report iterations
if     iter >=  0  && iter <  10
    fprintf(1,'    ---  iter =    %d;   abs res = %4.4e;   rel res = %4.4e \n',iter,resnorm,resnorm/resnorm0);
elseif iter >= 10  && iter < 100
    fprintf(1,'    ---  iter =   %d;   abs res = %4.4e;   rel res = %4.4e \n',iter,resnorm,resnorm/resnorm0);
elseif iter >= 100 && iter < 1000
    fprintf(1,'    ---  iter =  %d;   abs res = %4.4e;   rel res = %4.4e \n',iter,resnorm,resnorm/resnorm0);
end 
