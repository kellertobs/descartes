% record average phase velocities
HST.time(step) = time;

k = 0;
for it = 1:Nt
    HST.fp_sum  (step,it) = sum(fp(tp==it),  'all');     % phase fractions of sum of particles
    HST.fp_std  (step,it) = std(fp(tp==it),1,'all');     % phase fractions standard deviations (std)

    HST.Wp_mean(step,it) = mean(Wp(tp==it),  'all');     % particle settling z-speed mean
    HST.Wp_std (step,it) = std (Wp(tp==it),1,'all');     % particle settling z-speed std

    HST.Up_mean(step,it) = mean(Up(tp==it),  'all');     % particle settling x-speed mean
    HST.Up_std (step,it) = std (Up(tp==it),1,'all');     % particle settling x-speed std

    HST.Vp_mean(step,it) = mean(Vp(tp==it),  'all');     % particle settling speed mean

    HST.Wm_mean(step,it) = mean(Wm(tp==it),  'all');     % melt flow z-speed mean
    HST.Wm_std (step,it) = std (Wm(tp==it),1,'all');     % melt flow z-speed std

    HST.Um_mean(step,it) = mean(Um(tp==it),  'all');     % melt flow x-speed mean
    HST.Um_std (step,it) = std (Um(tp==it),1,'all');     % melt flow x-speed std

    HST.Vm_mean(step,it) = mean(Vm(tp==it),  'all');     % melt flow speed mean

    dWp = Wp - HST.Wp_mean(step,it);
    dUp = Up - 0;
    dVp = sqrt(dWp.^2 + dUp.^2);

    HST.dWp_mean(step,it) = mean(abs(dWp(tp==it)),'all');
    HST.dUp_mean(step,it) = mean(abs(dUp(tp==it)),'all');
    HST.dVp_mean(step,it) = mean(abs(dVp(tp==it)),'all');

    HST.DWp_mean(step,it) = mean(DWp(tp==it),  'all');   % particle-melt z-speed difference mean
    HST.DWp_std (step,it) = std (DWp(tp==it),1,'all');   % particle-melt z-speed difference std 

    HST.DVp_mean(step,it) = mean(DVp(tp==it),  'all');   % particle-melt speed difference mean

    HST.FHp_mean(step,it) = mean(DWp(tp==it)./DW0(it),  'all');   % hindering factor mean
    HST.FHp_std (step,it) = std (DWp(tp==it)./DW0(it),1,'all');   % hindering factor std

    HST.Wc_mean(step) = mean(abs(W),  'all');            % convection z-speed mean
    HST.Wc_std (step) = std (abs(W),1,'all');            % convection z-speed std

    HST.Uc_mean(step) = mean(abs(U),  'all');            % convection z-speed mean
    HST.Uc_std (step) = std (abs(U),1,'all');            % convection z-speed std

    HST.Vc_mean(step) = mean(abs(Vel),  'all');          % convection speed mean
    HST.Vc_std (step) = std (abs(Vel),1,'all');          % convection speed std

    % hindered settling diagnostics averaged over latter half of simulation
    HST.fp_tavg(step,it)  = mean(HST.fp_sum  (end-floor(step/2):end,it));
    HST.Wp_tavg(step,it)  = mean(HST.Wp_mean (end-floor(step/2):end,it));
    HST.DWp_tavg(step,it) = mean(HST.DWp_mean(end-floor(step/2):end,it));
    HST.FHp_tavg(step,it) = mean(HST.FHp_mean(end-floor(step/2):end,it));
    HST.Wc_tavg(step   )  = mean(HST.Wc_mean (end-floor(step/2):end   ));

    % particle positions
    for ip = 1:Np(it)
        k = k+1;
        HST.pos(step,k,1) = zp(k);
        HST.pos(step,k,2) = xp(k);
    end
end

