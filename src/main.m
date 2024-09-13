% initialise model run
init;

% physical time stepping loop
step = 1;
while time <= tend && step <= Nstep
    
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

    % write history
    history;

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
