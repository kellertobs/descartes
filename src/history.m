% record average phase velocities
HST.time(step) = time;

for it = 1:Nt
    HST.fp_sum  (step,it) = sum(fp(tp==it),  'all');     % phase fractions of sum of particles
    HST.fp_std  (step,it) = std(fp(tp==it),1,'all');     % phase fractions standard deviations (std)

    HST.Wp_mean(step,it) = mean(Wp(tp==it),  'all');     % particle velocities mean
    HST.Wp_std (step,it) = std (Wp(tp==it),1,'all');     % particle velocities std

    HST.DWp_mean(step,it) = mean(DWp(tp==it),  'all');   % particle segregation speeds mean
    HST.DWp_std (step,it) = std (DWp(tp==it),1,'all');   % particle segregation speeds std 

    HST.FHp_mean(step,it) = mean(DWp(tp==it)./DW0(it),  'all');   % hindering factor mean
    HST.FHp_std (step,it) = std (DWp(tp==it)./DW0(it),1,'all');   % hindering factor std

    % hindered settling diagnostics averaged over latter half of simulation
    HST.fp_tavg(step,it)  = mean(HST.fp_sum  (end-floor(step/2):end,it));
    HST.Wp_tavg(step,it)  = mean(HST.Wp_mean (end-floor(step/2):end,it));
    HST.DWp_tavg(step,it) = mean(HST.DWp_mean(end-floor(step/2):end,it));
    HST.FHp_tavg(step,it) = mean(HST.FHp_mean(end-floor(step/2):end,it));
end

