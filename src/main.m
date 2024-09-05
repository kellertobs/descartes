% initialise model run
init;

% physical time stepping loop
while time <= tend && step <= Nt
    
    % time step info
    timing;

    % store previous solution
    store;

    % solve fluid-mechanics equations
    fluidmech;

    % update non-linear parameters and auxiliary variables
    update;

    % report convergence
    report;

    % record model history
    history;

    % print model diagnostics
    diagnose;

    % plot model results
    if ~mod(step,nop); output; end
    
    % increment time/step
    time = time+dt;
    step = step+1;
    
end

% save final state of model
output;

diary off
