function parte()

    % time parameters
    dt = 0.001;
    sim_time = 1;
    t = 0:dt:sim_time;
    npts = length(t);
    
    % output storage
    dxs = zeros(1, npts);
    xs = zeros(1, npts);
    
    % simulation
    x = 0;
    dx = 0;
    x_f = 1;
    for i=1:npts
        dxs(i) = dx;
        xs(i) = x;
        [dx, x] = step(x_f, x, dx, dt);
    end
    
    % plots
    clf;
    subplot(1, 2, 1);    
    step_trace = ones(1, npts);
    step_trace(1) = 0;
    hold on;
    plot(t, step_trace, 'LineWidth', 1.5);
    plot(t, xs, 'LineWidth', 1.5);
    ylim([0 1.2]);
    xlim([-0.2 sim_time]);
    title('MainScope');
    xlabel('time');
    legend('x_f', 'x', 'Location', 'Southeast');
    
    subplot(1, 2, 2);
    plot(t, dxs, 'LineWidth', 1.5);
    ylim([0 5]);
    xlim([-0.2 sim_time]);
    title('DerivScope');
    xlabel('time');
    legend('dx', 'Location', 'Southeast');
end

function [dx_nplus1, x_nplus1] = step(x_f, x_n, dx_n, dt)
    alpha = 25;
    tau = 1;
    beta = 6;
    ddx = (alpha / tau) * (beta * (x_f - x_n) - dx_n);
    dx_nplus1 = dx_n + dt * ddx;
    x_nplus1 = x_n + dt * dx_nplus1;
end