function parte()

    % time parameters
    dt = 0.001;
    tau = 1;
    t = 0:dt:tau;
    npts = length(t);
    
    % output storage
    dxs = zeros(1, npts);
    xs = zeros(1, npts);
    
    % simulation
    x = 0;
    dx = 0;
    ddx = 0;
    x_f = 1;
    dx_f = 0;
    ddx_f = 0;
    for i=1:npts
        dxs(i) = dx;
        xs(i) = x;
        p = poly5(x, dx, ddx, x_f, dx_f, ddx_f, tau - t(i));
        dp = polyder(p);
        ddp = polyder(dp);
        x = polyval(p, dt);
        dx = polyval(dp, dt);
        ddx = polyval(ddp, dt);
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
    xlim([-0.2 tau]);
    title('MainScope');
    xlabel('time');
    legend('x_f', 'x', 'Location', 'Southeast');
    
    subplot(1, 2, 2);
    plot(t, dxs, 'LineWidth', 1.5);
    ylim([0 2.5]);
    xlim([-0.2 tau]);
    title('DerivScope');
    xlabel('time');
    legend('dx', 'Location', 'Northeast');
end

% outputs a vector of polynomial coefficients in descending order,
% suitable for Matlab's built-in polyval and polyder functions
function p = poly5(x0, dx0, ddx0, xf, dxf, ddxf, T)
	c3 = (-12*dx0*T - 8*dxf*T - 3*ddx0*T^2 + ddxf*T^2 - 20*x0 + 20*xf)/(2*T^3);
	c4 = (16*dx0*T + 14*dxf*T + 3*ddx0*T^2 - 2*ddxf*T^2 + 30*x0 - 30*xf)/(2*T^4);
	c5 = (-6*dx0*T - 6*dxf*T - ddx0*T^2 + ddxf*T^2 - 12*x0 + 12*xf)/(2*T^5);
	p  = [c5 c4 c3 ddx0/2 dx0 x0];
end
