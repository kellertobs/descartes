% record average phase velocities
HST.time(step) = time;

% Wpc = zeros(Nz/Ns,Nx/Ns,Nt);
% Fpc = zeros(Nz/Ns,Nx/Ns,Nt);
% Wmc = zeros(Nz/Ns,Nx/Ns,1 );
% Fmc = zeros(Nz/Ns,Nx/Ns,1 );
% 
% for i=1:Nx/Ns
%     for j=1:Nz/Ns
%         for it = 1:Nt
%             indx = (i-1)*Ns+(1:Ns);
%             indz = (j-1)*Ns+(1:Ns);
%             Wpc(j,i,it) = sum( Wc(indz,indx).*(C(indz,indx,it)>=0) ,'all')./sum( C(indz,indx,it)>=0 ,'all');
%             Fpc(j,i,it) = mean( C(indz,indx,it)>=0 ,'all');
%             if isnan(Wpc(j,i,it)); Wpc(j,i,it) = 0; end
%         end
%         Wmc(j,i) = sum( Wc(indz,indx).*all(C(indz,indx,:)<0,3) ,'all')./sum( all(C(indz,indx,:)<0,3) ,'all');
%         Fmc(j,i) = mean( all(C(indz,indx,:)<0,3) ,'all');
%     end
% end
% DWp = Wpc - Wmc;

for it = 1:Nt
    HST.DWp_NM (step,it) = mean( DWp(tp==it) ,  'all');
    HST.DWp_std(step,it) = std ( DWp(tp==it) ,1,'all');
end
HST.DWp_EM(step,:) = (rhop-mean(rho(:))).*grav.*rp.^2./geomean(eta(:));
HST.DWp_ST(step,:) = (rhop-rhom).*grav.*rp.^2./etam;
HST.DWp_HS(step,:) = (rhop-rhom).*grav.*rp.^2./etam.*(1-fp).^5;
HST.FHS(step,:)    = HST.DWp_NM(step,:)./HST.DWp_ST(step,:);