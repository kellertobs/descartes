% initialise model run
init;

% physical time stepping loop
step = 1;
while time <= tend && step <= Nt
    
    % time step info
    timing;

    % store previous solution
    store;

    % nonlinear iteration loop
    resnorm = 1; resnorm0 = 1;
    iter    = 1;
    while resnorm > atol && resnorm/resnorm0 > rtol && iter<=maxit

        % solve fluid-mechanics equations
        fluidmech;

        % update non-linear parameters and auxiliary variables
        update;

        % report convergence
        report;

        iter = iter+1;
    end

    % record average phase velocities
    HST.time(step) = time;
    Wc = -(W(1:end-1,2:end-1)+W(2:end,2:end-1))/2;
    for it = 1:Np
        HST.Wp(step,it) = mean(Wc(C(:,:,it)>=0));
    end
    HST.Wp(step,Np+1) = mean(Wc(any(C<0,3)));

    for it = 1:Np
        HST.stdWp(step,it) = std(Wc(C(:,:,it)>=0));
    end
    HST.Wp(step,Np+1) = std(Wc(any(C<0,3)));

    HST.DWp(step,:) = HST.Wp(step,1:Np) - HST.Wp(step,Np+1);
    HST.DWp0(step,:) = -(rhop-mean(rho(:))).*grav.*rp.^2./geomean(eta(:));

    % print model diagnostics
    % diagnose;

    % plot model results
    if ~mod(step,nop); output; end
    
    % increment time/step
    time = time+dt;
    step = step+1;
    if frst; frst = 0; end
end

% save final state of model
output;

diary off
