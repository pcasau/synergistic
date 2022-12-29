global z_0 r Gamma_1 theta_0 epsilon theta Gamma_2 k_u
z_0   = [1 0]';            % Center of the obstacle
r     = 0.5;               % Radius of the obstacle
Gamma_1 = eye(2);          % Estimator gain
Gamma_2 = eye(2);          % Backstepping gain
k_u = 1;                   % Derivative gain
theta_0 = 1;               % Maximum disturbance magnitude 
epsilon = 0.5;             % Proj operator parameter
theta = sqrt(2)*[1 1]'./2; % Disturbance

% Computation of partial derivatives
if ~exist('Dphi','file')
    partialDiffFnc(@(xi) phi(xi(1:3),xi(4)),[4,1],1:3,'Dphi',{1:3,4});
end
if ~exist('DV_0','file')
    partialDiffFnc(@(xi) V_0(xi(1:3),xi(4)),[4,1],1:3,'DV_0',{1:3,4});
end
if ~exist('Dxkappa_1','file')
    partialDiffFnc(@(xi) kappa_1(xi(1:3),xi(4:6)),[6,1],1:3,'Dxkappa_1',{1:3,4:6});
end
if ~exist('Dqkappa_1','file')
    partialDiffFnc(@(xi) kappa_1(xi(1:3),xi(4:6)),[6,1],4,'Dqkappa_1',{1:3,4:6});
end

% Simulation parameters
TSPAN = [0 10];            % Continuous-time bounds
JSPAN = [0 10];            % Discrete-time bounds
rule  = 1;                 % Priority given to jumps
z_i   = [2;0];             % Initial position on the plane
hat_theta_i = zeros(2,1);  % Initial theta estimate
u0 = zeros(2,1);           % Initial control input

%Simulation of the nominal system
%[t,j,xi] = HyEQsolver(@F_cl0, @G_cl0, @C_0, @D_0, xi_0, TSPAN, JSPAN, rule);

% Simulation of the perturbed system
hf = figure;
hf(2) = figure;
colors = get(0,'DefaultAxesColorOrder');
colors = colors([1 3],:);
linestyles = {'--','-'};
hp = [];
for J = 1:2
    for q_i = -1:2:1
        if J == 1
            xi0 = [psi(z_i);q_i;hat_theta_i];
        else
            V2_0 = V_2(psi(z_i),[q_i;hat_theta_i;u0],theta);
            xi0 = [psi(z_i);q_i;hat_theta_i;u0;V2_0];
        end
        eval(['[t,j,xi] = HyEQsolver(@F_cl' num2str(J) ...
            ', @G_Omega' num2str(J) ...
            ', @C_Omega' num2str(J) ...
            ', @D_Omega' num2str(J) ...
            ', xi0, TSPAN, JSPAN, rule);']);
        z = zeros(2,numel(t));
        for I = 1:numel(t)
            z(:,I) = psi_inv(xi(I,1:3)');
        end
        figure(hf(1))
        aux = plot(z(1,:),z(2,:),...
            'color',colors(J,:),'linestyle',linestyles{(q_i+3)/2});
        hp = [hp aux];
        hold on
        figure(hf(2))
        subplot(2,1,1)
        plot(t,sqrt(sum(z.^2,1)),...
            'color',colors(J,:),'linestyle',linestyles{(q_i+3)/2})
        hold on
        subplot(2,1,2)
        plot(t,sqrt(sum((xi(:,5:6)-repmat(theta',[numel(t),1])).^2,2)),...
            'color',colors(J,:),'linestyle',linestyles{(q_i+3)/2})
        hold on
    end
end

% Format plots
figure(hf(1))
grid on
xlabel('z_1')
ylabel('z_2')
axis equal
draw_circle('center',z_0,'radius',epsilon)
text(z0(1),z0(2),'N','horizontalalignment','center','verticalalignment','middle')
legend(hp,{'Section VI-B with q(0,0)=-1',...
    'Section VI-B with q(0,0)=1',...
    'Section VI-C with q(0,0)=-1',...
    'Section VI-C with q(0,0)=1'});
hold off

figure(hf(2))
subplot(2,1,1)
grid on
legend(hp,{'Section VI-B with q(0,0)=-1',...
    'Section VI-B with q(0,0)=1',...
    'Section VI-C with q(0,0)=-1',...
    'Section VI-C with q(0,0)=1'});
ylabel('|z(t)|')
subplot(2,1,2)
grid on
ylabel('Estimation error')
xlabel('t')


% Problem Setup
function x = psi(z)
    global z_0 r
    x = [log(norm(z-z_0)-r); (z-z_0)/norm(z-z_0)];
end

function z = psi_inv(x)
    global z_0 r
    z = z_0+x(2:3)*(exp(x(1))+r);
end

function out = Dpsi(z)
    global z_0 r
    out= [1/(norm(z-z_0)-r)*(z-z_0)'/norm(z-z_0);PT((z-z_0)/norm(z-z_0))/norm(z-z_0)];
end

function out = PT(x)
    n = numel(x);
    out = eye(n)-x*x';
end

function dx = F(x,q,u,theta)
    dx = Dpsi(psi_inv(x))*(u+theta);
end

% The Nominally Synergistic Controller
function out = U(x,q)
    if q*x(3) ~= 1
        out = 1;
    else
        out = 0;
    end
end

function dq = F_c(x,q)
    dq = 0;
end

function coords = phi(x,q)
    coords = [x(1),x(2)/(1-q*x(3))]';
end

function out = V_0(x,q)
    if U(x,q)
        out = 0.5*norm(phi(x,q)-phi(psi(zeros(2,1)),q))^2;
    else
        out = Inf;
    end
end

function out = kappa_0(x,q)
    out = -Dpsi(psi_inv(x))'*Dphi(x,q)'*(phi(x,q)-phi(psi(zeros(2,1)),q));
end

function out = delta(x,q)
    out = 1;
end

function out = nu_V0(x,q)
    out = min([V_0(x,q),V_0(x,-q)]);
end

function out = varrho_V0(x,q)
    if V_0(x,q)-V_0(x,-q) > 0
        out = -q;
    elseif V_0(x,q)-V_0(x,-q) < 0
        out = q;
    else
        % The dynamics are set-valued
        %out = datasample([-1 1]);  
        % For reproducibility purposes we default to:
        out = q;
    end
end

function out = mu_V0(x,q)
    out = V_0(x,q)-nu_V0(x,q);
end

function dxi = F_cl0(xi)
    x = xi(1:3);
    q = xi(4);
    theta = zeros(2,1); % Since this is the nominal case
    dx = F(x,q,kappa_0(x,q),theta);
    dq = F_c(x,q);
    dxi = [dx;dq];
end

function out = C_0(xi)
    x = xi(1:3);
    q = xi(4);
    if mu_V0(x,q) <= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function out = D_0(xi)
    x = xi(1:3);
    q = xi(4);
    if mu_V0(x,q) >= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function next_xi = G_cl0(xi)
    x = xi(1:3);
    q = xi(4);
    next_xi = [x;varrho_V0(x,q)];
end

% Adaptive Backstepping of Synergistic Hybrid Feedback for Affine Control Systems
function out = Proj(eta, hat_theta)
    global  theta_0 epsilon
    ell = numel(eta);
    p = (hat_theta'*hat_theta-theta_0^2)/(epsilon^2+2*epsilon*theta_0);
    gradp = 2*hat_theta/(epsilon^2+2*epsilon*theta_0);
    if p<=0 || gradp'*eta<=0
        out = eta;
    else
        out = (eye(ell)-p*(gradp*gradp')/norm(gradp)^2)*eta;
    end
end

function u = kappa_1(x,x_c1)
    q = x_c1(1);
    hat_theta = x_c1(2:3);
    u = kappa_0(x,q)-hat_theta;
end

function out = V_1(x,x_c1,theta)
    global Gamma_1
    q = x_c1(1);
    hat_theta = x_c1(2:3);
    out = V_0(x,q)+0.5*norm(Gamma_1^(-1/2)*(theta-hat_theta))^2;
end

function dx_c1 = F_c1(x,x_c1)
    global Gamma_1 
    q = x_c1(1);
    hat_theta = x_c1(2:3);
    dq = F_c(x,q);
    dhat_theta = Gamma_1*Proj(Dpsi(psi_inv(x))'*DV_0(x,q)',hat_theta);
    dx_c1 = [dq;dhat_theta];
end

function out = nu_V1(x,x_c1,theta)
    q = x_c1(1);
    out = nu_V0(x,q);
end

function out = varrho_V1(x,x_c1,theta)
    q = x_c1(1);
    hat_theta = x_c1(2:3);
    out = [varrho_V0(x,q);
        theta];
end

function out = mu_V1(x,x_c1,theta)
    global Gamma_1
    q = x_c1(1);
    hat_theta = x_c1(2:3);
    out = mu_V0(x,q)+0.5*norm(Gamma_1^{-1/2}*(theta-hat_theta))^2;
end

function out = min_mu_V1(x,x_c1)
    global theta_0
    q = x_c1(1);
    hat_theta = x_c1(2:3);
    if norm(hat_theta) <= theta_0
        out = mu_V0(x,q);
    else
        % These computations assume Gamma_1 = eye(2)
        out = mu_V0(x,q)+0.5*(1-theta_0/norm(hat_theta))^2*norm(hat_theta)^2;
    end
end

function out = hat_G(hat_theta)
    global theta_0
    if norm(hat_theta) <= theta_0
        out = hat_theta;
    else
        out = theta_0*hat_theta/norm(hat_theta);
    end
end

function out = G_c1(x,x_c1)
    q = x_c1(1);
    hat_theta = x_c1(2:3);
    out = [varrho_V0(x,q);
        hat_G(hat_theta)];
end

function out = C_Omega1(xi)
    x = xi(1:3);
    x_c1 = xi(4:6);
    q = x_c1(1);
    if min_mu_V1(x,x_c1) <= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function dxi = F_cl1(xi)
    global theta
    x = xi(1:3);
    x_c1 = xi(4:6);
    q = x_c1(1);
    dx = F(x,q,kappa_1(x,x_c1),theta);
    dx_c1 = F_c1(x,x_c1);
    dxi = [dx;dx_c1];
end

function out = D_Omega1(xi)
    x = xi(1:3);
    x_c1 = xi(4:6);
    q = x_c1(1);
    if min_mu_V1(x,x_c1) >= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function next_xi = G_Omega1(xi)
    x = xi(1:3);
    x_c1 = xi(4:6);
    next_x = x;
    next_x_c1 = G_c1(x,x_c1);
    next_xi = [next_x;
        next_x_c1];
end

% Backstepping
function out = v(x,x_c2)
    global Gamma_2
    x_c1 = x_c2(1:3);
    q = x_c1(1);
    u = x_c2(4:5);
    out = Dpsi(psi_inv(x))'*(DV_0(x,q)'-Dxkappa_1(x,x_c1)'*Gamma_2^(-1)*(u-kappa_1(x,x_c1)));
end

function out = f_u(x,x_c2)
    global Gamma_1 Gamma_2 k_u
    x_c1 = x_c2(1:3);
    q = x_c1(1);
    hat_theta = x_c1(2:3);
    u = x_c2(4:5);
    out = -Gamma_1*Proj(v(x,x_c2),hat_theta)-k_u*(u-kappa_1(x,x_c1))...
        -Gamma_2*Dpsi(psi_inv(x))'*DV_0(x,q)'+Dxkappa_1(x,x_c1)*F(x,q,u,hat_theta);
end

function dx_c2 = F_c2(x,x_c2)
    global Gamma_1
    x_c1 = x_c2(1:3);
    q = x_c1(1);
    hat_theta = x_c1(2:3);
    u = x_c2(4:5);
    dx_c2 = [F_c(x,q);
        Gamma_1*Proj(v(x,x_c2),hat_theta);
        f_u(x,x_c2)+Dqkappa_1(x,x_c1)*F_c(x,q)];
end

function out = V_2(x,x_c2,theta)
    global Gamma_2
    x_c1 = x_c2(1:3);
    u = x_c2(4:5);
    out = V_1(x,x_c1,theta)+0.5*norm(Gamma_2^(-1/2)*(u-kappa_1(x,x_c1)))^2;
end

function g_c2 = G_c2(x,x_c2)
    x_c1 = x_c2(1:3);
    g_c1 = G_c1(x,x_c1);
    g_u = kappa_1(x,x_c1);
    g_c2 = [g_c1;g_u];
end

function out = min_mu_V2(x,x_c2)
    global Gamma_2
    x_c1 = x_c2(1:3);
    u = x_c2(4:5);
    out = min_mu_V1(x,x_c1)+0.5*norm(Gamma_2^(-1/2)*(u-kappa_1(x,x_c1)))^2;
end

function out = C_Omega2(xi)
    x = xi(1:3);
    x_c2 = xi(4:8);
    q = x_c2(1);
    if min_mu_V2(x,x_c2) <= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function dxi = F_cl2(xi)
    global theta Gamma_2 k_u
    x = xi(1:3);
    x_c2 = xi(4:8);
    x_c1 = x_c2(1:3);
    q = x_c2(1);
    hat_theta = x_c2(2:3);
    u = x_c2(4:5);
    dx = F(x,q,u,theta);
    dx_c2 = F_c2(x,x_c2);
    dV2 = DV_0(x,q)*F(x,q,kappa_0(x,q),zeros(2,1))...
        -k_u*(u-kappa_1(x,x_c1))'*Gamma_2^(-1)*(u-kappa_1(x,x_c1));
    dxi = [dx;dx_c2;dV2];
end

function out = D_Omega2(xi)
    x = xi(1:3);
    x_c2 = xi(4:8);
    q = x_c2(1);
    if min_mu_V2(x,x_c2) >= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function next_xi = G_Omega2(xi)
    global theta
    x = xi(1:3);
    x_c2 = xi(4:8);
    next_x = x;
    next_x_c2 = G_c2(x,x_c2);
    next_xi = [next_x;next_x_c2;V_2(next_x,next_x_c2,theta)];
end