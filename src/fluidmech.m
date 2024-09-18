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

% Left boundary (periodic)
ii  = MapW(:,1); jj1 = ii; jj2 = MapW(:,end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% Right boundary (periodic)
ii  = MapW(:,end); jj1 = ii; jj2 = MapW(:,2);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% Internal points
ii    = MapW(:,2:end-1);
EtaC1 = etaco(:     ,1:end-1);  EtaC2 = etaco(:     ,2:end);
EtaP1 = etacc(icz(1:end-1),:);  EtaP2 = etacc(icz(2:end),:);

% Coefficients multiplying z-velocities W
%             top              ||          bottom             
jj1 = MapW(ifz(1:end-2),2:end-1); jj2 = MapW(ifz(3:end),2:end-1); 
%             left             ||          right
jj3 = MapW(ifz(2:end-1),1:end-2); jj4 = MapW(ifz(2:end-1),3:end);

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
jj1 = MapU(1:end-1,1:end-1); jj2 = MapU(2:end,1:end-1); 
%        top right        ||       bottom right
jj3 = MapU(1:end-1,2:end  ); jj4 = MapU(2:end,2:end  );

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % U one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % U one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and right


% -------------------------------------------------------------------------
% Assemble coefficients of X-Stress Divergence (X-momentum)
% -------------------------------------------------------------------------

% Top boundary (periodic)
ii  = MapU(1,:); jj1 = ii; jj2 = MapU(end-1,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)]; AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)]; AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% Bottom boundary (periodic)
ii  = MapU(end,:); jj1 = ii; jj2 = MapU(2,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)]; AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)]; AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% Internal points
ii    = MapU(2:end-1,:);
EtaC1 = etaco(1:end-1,:     );  EtaC2 = etaco(2:end,:     );
EtaP1 = etacc(:,icx(1:end-1));  EtaP2 = etacc(:,icx(2:end));

% Coefficients multiplying x-velocities U
%            left              ||          right         
jj1 = MapU(2:end-1,ifx(1:end-2)); jj2 = MapU(2:end-1,ifx(3:end)); 
%            top               ||          bottom
jj3 = MapU(1:end-2,ifx(2:end-1)); jj4 = MapU(3:end,ifx(2:end-1));

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
jj1 = MapW(1:end-1,1:end-1); jj2 = MapW(1:end-1,2:end); 
%       bottom left       ||       bottom right
jj3 = MapW(2:end  ,1:end-1); jj4 = MapW(2:end  ,2:end);

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
ii  = MapW(:,2:end-1);
advn_mz = advect(rhofz.*W(:,2:end-1),(U(1:end-1,:)+U(2:end,:))/2,(W(ifz(1:end-1),2:end-1)+W(ifz(2:end),2:end-1))/2,h,{ADVN,''},[1,2],BCA);
advn_mz([1 end],:) = repmat((advn_mz(1,:)+advn_mz(end,:))/2,2,1);
rr  = + (rhofz - mean(rhofz,2)) .* grav ...
      + (a2.*rhoWo+a3.*rhoWoo)/dt ...
      - advn_mz;

IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];

% X-RHS vector
ii  = MapU(2:end-1,:);
advn_mx = advect(rhofx.*U(2:end-1,:),(U(2:end-1,ifx(1:end-1))+U(2:end-1,ifx(2:end)))/2,(W(:,1:end-1)+W(:,2:end))/2,h,{ADVN,''},[1,2],BCA);
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
    ii  = MapW(:,2:end-1);
    jj1 = MapP(1:end-1,2:end-1); jj2 = MapP(2:end,2:end-1);

    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)]; AAL = [AAL; aa(:)-1/h];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)]; AAL = [AAL; aa(:)+1/h];

    % X-gradient
    ii  = MapU(2:end-1,:);
    jj1 = MapP(2:end-1,1:end-1); jj2 = MapP(2:end-1,2:end);

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
    ii = MapP(2:end-1,2:end-1);

    % U, W velocities
    jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end);
    jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);

    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)]; AAL = [AAL; aa(:)-1/h];  % U left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)]; AAL = [AAL; aa(:)+1/h];   % U right
    IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)]; AAL = [AAL; aa(:)-1/h];  % W top
    IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)]; AAL = [AAL; aa(:)+1/h];   % W bottom

    % Assemble DD matrix
    DD = sparse(IIL, JJL, AAL, NP, NW + NU);

end


% -------------------------------------------------------------------------
% Assemble Pressure Matrix and RHS (KP, RP)
% -------------------------------------------------------------------------

if ~exist('KP', 'var')

    IIL = []; JJL = []; AAL = [];  % Reset for KP assembly

    % boundary points
    ii  = [MapP(1    ,:); MapP(end,:)]; % top & bot
    jj1 = ii;
    jj2 = [MapP(end-1,:); MapP(2  ,:)];

    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];

    ii  = [MapP(:,1    ); MapP(:,end)]; % left & right
    jj1 = ii;
    jj2 = [MapP(:,end-1); MapP(:,2  )];

    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];

    % Assemble KP matrix
    KP = sparse(IIL, JJL, AAL, NP, NP);

end

% RHS assembly for pressure (RP)
IIR = []; AAR = [];
ii  = MapP(2:end-1,2:end-1);
rr  = zeros(size(ii));
IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];

RP = sparse(IIR, ones(size(IIR)), AAR, NP, 1);

% Set P = 0 at fixed point (pressure anchor)
nzp = round((Nz+2)/2);  % Midpoint in Z
nxp = round((Nx+2)/2);  % Midpoint in X
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
W = full(reshape(SOL(MapW(:))        ,Nz+1,Nx+2));  % matrix z-velocity
U = full(reshape(SOL(MapU(:))        ,Nz+2,Nx+1));  % matrix x-velocity
P = full(reshape(SOL(MapP(:)+(NW+NU)),Nz+2,Nx+2));  % matrix dynamic pressure
P = P - mean(P(:));

% Calculate velocity magnitude for diagnostics
Wc  = (W(1:end-1,2:end-1) + W(2:end,2:end-1)) / 2;
Uc  = (U(2:end-1,1:end-1) + U(2:end-1,2:end)) / 2;
Vel = sqrt(Wc.^2 + Uc.^2);




% % assemble coefficients for matrix velocity diagonal and right-hand side
% 
% IIL = [];       % equation indeces into L
% JJL = [];       % variable indeces into L
% AAL = [];       % coefficients for L
% IIR = [];       % equation indeces into R
% AAR = [];       % forcing entries for R


% % assemble coefficients of z-stress divergence
% 
% % left boundary
% ii  = MapW(:,1); jj1 = ii; jj2 = MapW(:,end-1);
% aa  = zeros(size(ii));
% IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
% IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
% IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];
% 
% % right boundary
% ii  = MapW(:,end); jj1 = ii; jj2 = MapW(:,2);
% aa  = zeros(size(ii));
% IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
% IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
% IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];
% 
% % internal points
% ii    = MapW(:,2:end-1);
% EtaC1 = etaco(:     ,1:end-1);  EtaC2 = etaco(:     ,2:end);
% EtaP1 = etacc(icz(1:end-1),:);  EtaP2 = etacc(icz(2:end),:);
% 
% % coefficients multiplying z-velocities W
% %             top          ||         bottom          ||           left            ||          right
% jj1 = MapW(ifz(1:end-2),2:end-1); jj2 = MapW(ifz(3:end),2:end-1); jj3 = MapW(ifz(2:end-1),1:end-2); jj4 = MapW(ifz(2:end-1),3:end);
% 
% aa  = a1.*rhofz./dt;
% IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term
% 
% aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
% IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % W on stencil centre
% IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % W one above
% IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % W one below
% IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % W one to the left
% IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % W one to the right
% 
% % coefficients multiplying x-velocities U
% %         top left         ||        bottom left          ||       top right       ||       bottom right
% jj1 = MapU(1:end-1,1:end-1); jj2 = MapU(2:end,1:end-1); jj3 = MapU(1:end-1,2:end); jj4 = MapU(2:end,2:end);
% 
% IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % U one to the top and left
% IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and left
% IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % U one to the top and right
% IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and right
% 
% 
% % z-RHS vector
% advn_mz = advect(rhofz.*W(:,2:end-1),(U(1:end-1,:)+U(2:end,:))/2,(W(ifz(1:end-1),2:end-1)+W(ifz(2:end),2:end-1))/2,h,{ADVN,''},[1,2],BCA);
% advn_mz([1 end],:) = repmat((advn_mz(1,:)+advn_mz(end,:))/2,2,1);
% rr  = + (rhofz - mean(rhofz,2)) .* grav ...
%       + (a2.*rhoWo+a3.*rhoWoo)/dt ...
%       - advn_mz;
% 
% IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];
% 
% 
% %  assemble coefficients of x-stress divergence

% % top boundary
% ii  = MapU(1,:); jj1 = ii; jj2 = MapU(end-1,:);
% aa  = zeros(size(ii));
% IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
% IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
% IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];
% 
% % bottom boundary
% ii  = MapU(end,:); jj1 = ii; jj2 = MapU(2,:);
% aa  = zeros(size(ii));
% IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
% IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
% IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];
% 
% % internal points
% ii    = MapU(2:end-1,:);
% EtaC1 = etaco(1:end-1,:     );  EtaC2 = etaco(2:end,:     );
% EtaP1 = etacc(:,icx(1:end-1));  EtaP2 = etacc(:,icx(2:end));

% % coefficients multiplying x-velocities U
% %            left          ||          right          ||           top             ||          bottom
% jj1 = MapU(2:end-1,ifx(1:end-2)); jj2 = MapU(2:end-1,ifx(3:end)); jj3 = MapU(1:end-2,ifx(2:end-1)); jj4 = MapU(3:end,ifx(2:end-1));
% 
% aa  = a1.*rhofx./dt;
% IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term
% 
% aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
% IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % U on stencil centre
% IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % U one to the left
% IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % U one to the right
% IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % U one above
% IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % U one below
% 
% 
% % coefficients multiplying z-velocities W
% %         top left         ||        top right          ||       bottom left       ||       bottom right
% jj1 = MapW(1:end-1,1:end-1); jj2 = MapW(1:end-1,2:end); jj3 = MapW(2:end,1:end-1); jj4 = MapW(2:end,2:end);
% 
% IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
% IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the top and right
% IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the bottom and left
% IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right
% 
% % x-RHS vector
% advn_mx = advect(rhofx.*U(2:end-1,:),(U(2:end-1,ifx(1:end-1))+U(2:end-1,ifx(2:end)))/2,(W(:,1:end-1)+W(:,2:end))/2,h,{ADVN,''},[1,2],BCA);
% advn_mx(:,[1 end]) = repmat((advn_mx(:,1)+advn_mx(:,end))/2,1,2);
% rr  = + (a2.*rhoUo+a3.*rhoUoo)/dt ...
%       - advn_mx;
% 
% IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];
% 
% 
% % assemble coefficient matrix & right-hand side vector
% KV  = sparse(IIL,JJL,AAL,NW+NU,NW+NU);
% RV  = sparse(IIR,ones(size(IIR)),AAR);
% 
% 
% % assemble coefficients for gradient operator
% 
% if ~exist('GG','var')
%     IIL = [];       % equation indeces into A
%     JJL = [];       % variable indeces into A
%     AAL = [];       % coefficients for A
% 
%     % coefficients for z-gradient
%     ii  = MapW(:,2:end-1);
% 
%     %         top              ||          bottom
%     jj1 = MapP(1:end-1,2:end-1); jj2 = MapP(2:end,2:end-1);
% 
%     aa  = zeros(size(ii));
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the top
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the bottom
% 
% 
%     % coefficients for x-gradient
%     ii  = MapU(2:end-1,:);
% 
%     %         left             ||           right
%     jj1 = MapP(2:end-1,1:end-1); jj2 = MapP(2:end-1,2:end);
%     aa  = zeros(size(ii));
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the left
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the right
% 
%     % assemble coefficient matrix
%     GG  = sparse(IIL,JJL,AAL,NW+NU,NP);
% end
% 
% 
% % assemble coefficients for divergence operator
% 
% if ~exist('DD','var')
%     IIL = [];       % equation indeces into A
%     JJL = [];       % variable indeces into A
%     AAL = [];       % coefficients for A
% 
%     %internal points
%     ii  = MapP(2:end-1,2:end-1);
% 
%     % coefficients multiplying velocities U, W
%     %          left U          ||           right U       ||           top W           ||          bottom W
%     jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);
% 
%     aa  = zeros(size(ii));
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];  % U one to the left
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];  % U one to the right
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; aa(:)-1/h];  % W one above
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; aa(:)+1/h];  % W one below
% 
%     % assemble coefficient matrix
%     DD  = sparse(IIL,JJL,AAL,NP,NW+NU);
% end
% 
% 
% % assemble coefficients for matrix pressure diagonal and right-hand side
% 
% if ~exist('KP','var')
%     IIL = [];       % equation indeces into A
%     JJL = [];       % variable indeces into A
%     AAL = [];       % coefficients for A
% 
%     % boundary points
%     ii  = [MapP(1    ,:); MapP(end,:)]; % top & bot
%     jj1 = ii;
%     jj2 = [MapP(end-1,:); MapP(2  ,:)];
% 
%     aa  = zeros(size(ii));
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
% 
%     ii  = [MapP(:,1    ); MapP(:,end)]; % left & right
%     jj1 = ii;
%     jj2 = [MapP(:,end-1); MapP(:,2  )];
% 
%     aa  = zeros(size(ii));
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
% 
%     % assemble coefficient matrix
%     KP  = sparse(IIL,JJL,AAL,NP,NP);
% end
% 
% % KV(NW/2,:) = KV(NW/2,:) + 1;
% 
% % RHS
% IIR = [];       % equation indeces into R
% AAR = [];       % forcing entries for R
% 
% ii  = MapP(2:end-1,2:end-1);
% 
% rr  = zeros(size(ii));
% 
% IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];
% 
% % assemble right-hand side vector
% RP  = sparse(IIR,ones(size(IIR)),AAR,NP,1);
% 
% % KP(end,:) = 1;
% 
% % set P = 0 in fixed point
% nzp = round((Nz+2)/2);
% nxp = round((Nx+2)/2);
% np0 = MapP(nzp,nxp);
% DD(np0,:  ) = 0;
% KP(np0,:  ) = 0;
% KP(np0,np0) = sqrt(h^2./geomean(eta(:)));
% RP(np0,:  ) = 0;
% 
% 
% % assemble and scale global coefficient matrix and right-hand side vector
% 
% LL  = [ KV   GG  ; ...
%         DD   KP ];
% 
% RR  = [RV; RP];
% 
% SCL = (abs(diag(LL))).^0.5;
% SCL = diag(sparse( 1./(SCL + sqrt(h^2./geomean(eta(:)))) ));
% 
% FF  = LL*[W(:);U(:);P(:)] - RR;
% 
% LL  = SCL*LL*SCL;
% FF  = SCL*FF;
% RR  = SCL*RR;
% 
% 
% % Solve linear system of equations for vx, vz, P
% 
% SOL = SCL*(LL\RR);  % update solution
% % SOL = (LL\RR);  % update solution
% 
% % map solution vector to 2D arrays
% W = full(reshape(SOL(MapW(:))        ,Nz+1,Nx+2));  % matrix z-velocity
% U = full(reshape(SOL(MapU(:))        ,Nz+2,Nx+1));  % matrix x-velocity
% P = full(reshape(SOL(MapP(:)+(NW+NU)),Nz+2,Nx+2));  % matrix dynamic pressure
% P = P - mean(P(:));
% 
% Wc  = (W(1:end-1,2:end-1)+W(2:end,2:end-1))/2;
% Uc  = (U(2:end-1,1:end-1)+U(2:end-1,2:end))/2;
% Vel = sqrt(Wc.^2 + Uc.^2);