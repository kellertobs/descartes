% FLUIDMECH - Main numerical solver for fluid mechanics in Descartes
% Assemble coefficients for matrix velocity and pressure diagonal, gradient, 
% and divergence block matrices, scale and solve the resulting linear system.


% -------------------------------------------------------------------------
% Initialize coefficient arrays for the matrix and right-hand side (RHS)
% -------------------------------------------------------------------------
IIL = [];  % Equation indices into L (global coefficient matrix)
JJL = [];  % Variable indices into L
AAL = [];  % Coefficients for L
IIR = [];  % Equation indices into RHS vector R
AAR = [];  % Forcing entries for R


% -------------------------------------------------------------------------
% Assemble coefficients of Z-Stress Divergence (Z-momentum)
% -------------------------------------------------------------------------

% Internal points
ii    = MapW;
EtaC1 = etaco(:     ,1:end-1);  EtaC2 = etaco(:     ,2:end);
EtaP1 = etacc(icz(1:end-1),:);  EtaP2 = etacc(icz(2:end),:);

% Coefficients multiplying z-velocities W
%             top              ||          bottom             
jj1 = MapW(ifz(1:end-2),icx(2:end-1)); jj2 = MapW(ifz(3:end),icx(2:end-1)); 
%             left             ||          right
jj3 = MapW(ifz(2:end-1),icx(1:end-2)); jj4 = MapW(ifz(2:end-1),icx(3:end));

% Inertial term
aa  = a1.*rhofz./dt;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term

% Z-stress divergence terms (central stencil + neighbors)
aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % W on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % W one below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % W one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % W one to the right

% Coefficients multiplying x-velocities U
%        top left         ||        bottom left       
jj1 = MapU(icz(1:end-1),1:end-1); jj2 = MapU(icz(2:end),1:end-1); 
%        top right        ||       bottom right
jj3 = MapU(icz(1:end-1),2:end  ); jj4 = MapU(icz(2:end),2:end  );

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % U one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % U one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and right


% -------------------------------------------------------------------------
% Assemble coefficients of X-Stress Divergence (X-momentum)
% -------------------------------------------------------------------------

% Internal points
ii    = MapU;
EtaC1 = etaco(1:end-1,:     );  EtaC2 = etaco(2:end,:     );
EtaP1 = etacc(:,icx(1:end-1));  EtaP2 = etacc(:,icx(2:end));

% Coefficients multiplying x-velocities U
%            left              ||          right         
jj1 = MapU(icz(2:end-1),ifx(1:end-2)); jj2 = MapU(icz(2:end-1),ifx(3:end)); 
%            top               ||          bottom
jj3 = MapU(icz(1:end-2),ifx(2:end-1)); jj4 = MapU(icz(3:end),ifx(2:end-1));

% Inertial term
aa  = a1.*rhofx./dt;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term

% X-stress divergence terms (central stencil + neighbors)
aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % U on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % U one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % U one below

% Coefficients multiplying z-velocities W
%       top left          ||        top right         
jj1 = MapW(1:end-1,icx(1:end-1)); jj2 = MapW(1:end-1,icx(2:end)); 
%       bottom left       ||       bottom right
jj3 = MapW(2:end  ,icx(1:end-1)); jj4 = MapW(2:end  ,icx(2:end));

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right

% Assemble coefficient matrix for velocity
KV = sparse(IIL, JJL, AAL, NW + NU, NW + NU);


% -------------------------------------------------------------------------
% Assemble Z- and X-RHS Vectors
% -------------------------------------------------------------------------

% Z-RHS vector
ii  = MapW;
advn_mz = advect(rhoW,(U(icz(1:end-1),:)+U(icz(2:end),:))/2,(W(ifz(1:end-1),icx(2:end-1))+W(ifz(2:end),icx(2:end-1)))/2,h,{ADVN,''},[1,2],BCA);
advn_mz([1 end],:) = repmat((advn_mz(1,:)+advn_mz(end,:))/2,2,1);
rr  = + (rhofz - mean(rhofz,2)) .* grav ...
      + (a2.*rhoWo+a3.*rhoWoo)/dt ...
      - advn_mz;

IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];

% X-RHS vector
ii  = MapU;
advn_mx = advect(rhoU,(U(icz(2:end-1),ifx(1:end-1))+U(icz(2:end-1),ifx(2:end)))/2,(W(:,icx(1:end-1))+W(:,icx(2:end)))/2,h,{ADVN,''},[1,2],BCA);
advn_mx(:,[1 end]) = repmat((advn_mx(:,1)+advn_mx(:,end))/2,1,2);
rr  = + (a2.*rhoUo+a3.*rhoUoo)/dt ...
      - advn_mx;

IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];

% Assemble RHS vector for velocity
RV = sparse(IIR, ones(size(IIR)), AAR);


% -------------------------------------------------------------------------
% Assemble Gradient Operator Coefficients (GG)
% -------------------------------------------------------------------------

if ~exist('GG', 'var')

    IIL = []; JJL = []; AAL = [];  % Reset for GG assembly

    % Z-gradient
    ii  = MapW;
    jj1 = MapP(icz(1:end-1),icx(2:end-1)); jj2 = MapP(icz(2:end),icx(2:end-1));

    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)]; AAL = [AAL; aa(:)-1/h];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)]; AAL = [AAL; aa(:)+1/h];

    % X-gradient
    ii  = MapU;
    jj1 = MapP(icz(2:end-1),icx(1:end-1)); jj2 = MapP(icz(2:end-1),icx(2:end));

    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)]; AAL = [AAL; aa(:)-1/h];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)]; AAL = [AAL; aa(:)+1/h];

    % Assemble GG matrix
    GG = sparse(IIL, JJL, AAL, NW + NU, NP);

end


% -------------------------------------------------------------------------
% Assemble Divergence Operator Coefficients (DD)
% -------------------------------------------------------------------------

if ~exist('DD', 'var')

    IIL = []; JJL = []; AAL = [];  % Reset for DD assembly

    % Internal points for divergence operator
    ii = MapP;

    % U, W velocities
    jj1 = MapU(icz(2:end-1),ifx(2:end-2)); jj2 = MapU(icz(2:end-1),ifx(3:end-1));
    jj3 = MapW(ifz(2:end-2),icx(2:end-1)); jj4 = MapW(ifz(3:end-1),icx(2:end-1));

    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)]; AAL = [AAL; aa(:)-1/h];  % U left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)]; AAL = [AAL; aa(:)+1/h];   % U right
    IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)]; AAL = [AAL; aa(:)-1/h];  % W top
    IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)]; AAL = [AAL; aa(:)+1/h];   % W bottom

    % Assemble DD matrix
    DD = sparse(IIL, JJL, AAL, NP, NW + NU);

end

KP = sparse(NP, NP);
RP = sparse(NP,  1);

% Set P = 0 at fixed point (pressure anchor)
nzp = round(Nz/2);  % Midpoint in Z
nxp = round(Nx/2);  % Midpoint in X
np0 = MapP(nzp,nxp);    % Fixed point in pressure
DD(np0,:)  = 0;
KP(np0,:)  = 0;
KP(np0,np0) = sqrt(h^2 / geomean(eta(:)));
RP(np0,:)  = 0;


% -------------------------------------------------------------------------
% Assemble, Scale, and Solve Global System
% -------------------------------------------------------------------------

% Assemble global matrix (LL) and RHS vector (RR)
LL = [KV GG; ...
      DD KP];

RR = [RV; RP];

% Scale system for numerical stability
SCL = sqrt(abs(diag(LL)));
SCL = diag(sparse(1 ./ (SCL + sqrt(h^2 / geomean(eta(:))))));

LL  = SCL*LL*SCL;
RR  = SCL*RR;

% Solve linear system
SOL = SCL * (LL \ RR);

% Map solution vector back to velocity and pressure fields
W = full(reshape(SOL(MapW(:))        ,Nz+1,Nx));  % matrix z-velocity
U = full(reshape(SOL(MapU(:))        ,Nz,Nx+1));  % matrix x-velocity
P = full(reshape(SOL(MapP(:)+(NW+NU)),Nz,Nx));  % matrix dynamic pressure
P = P - mean(P(:));

% Calculate velocity magnitude for diagnostics
Wc  = (W(1:end-1,:) + W(2:end,:)) / 2;
Uc  = (U(:,1:end-1) + U(:,2:end)) / 2;
Vel = sqrt(Wc.^2 + Uc.^2);

resnorm_FM = 0;