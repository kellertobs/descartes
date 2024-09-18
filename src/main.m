% Main routine for numerical modelling code Descartes created to simulate
% particle settling/flotation in a small periodic domain

% Initialize the model (this sets up initial conditions, parameters, etc.)
init;

% Physical time stepping loop
step = 1;
while time <= tend && step <= Nstep

    % Display timing information for the current step
    timing;
    
    % Store previous solution
    store;
    
    % Initialize nonlinear iteration loop variables
    resnorm = 1;         % Residual norm (for convergence check)
    resnorm0 = 1;        % Initial residual norm (for relative convergence)
    iter = 1;            % Nonlinear iteration counter
    
    % Nonlinear iteration loop (solve until convergence or max iterations)
    while resnorm > atol && resnorm/resnorm0 > rtol && iter <= maxit
        
        % Solve fluid-mechanics equations (fluidmech should be vectorized internally)
        fluidmech;
        
        % Update nonlinear parameters and auxiliary variables
        update;
        
        % Report convergence progress
        report;
        
        % Increment the nonlinear iteration counter
        iter = iter + 1;

    end
    
    % Record history of the current step
    history;
    
    % Optionally output model results at a fixed interval (every 'nop' steps)
    if ~mod(step, nop)
        output;
    end
    
    % Increment time and step
    time = time + dt;   % Advance time by the current time step
    step = step + 1;    % Increment step counter
    
    % Reset first-step after first step
    if frst
        frst = 0;
    end
end

% Save the final state of the model after the last step
output;

% Close logging/diary file if it was enabled
diary off;