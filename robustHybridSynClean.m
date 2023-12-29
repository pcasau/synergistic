global z_0 r Gamma_0 Gamma_1 theta_0 epsilon theta Gamma_2 k_u gamma hat_W
z_0   = [1 0]';            % Center of the obstacle
r     = 0.5;               % Radius of the obstacle
Gamma_0 = 2*eye(2);        % Nominal feedback gain
Gamma_1 = eye(2);          % Estimator gain
Gamma_2 = eye(2);          % Backstepping gain
k_u = 4;                   % Derivative gain
theta_0 = 1;               % Maximum disturbance magnitude 
epsilon = 0.5;             % Proj operator parameter
theta = sqrt(2)*[1 1]'./4; % Disturbance
gamma = 0.5;               % X_c: q'*x0 <= gamma 
hat_W = eye(2);            % Matched uncertainty - See Assumption 2

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
hf(3) = figure;
colors = get(0,'DefaultAxesColorOrder');
colors = colors([1 3],:);
linestyles = {'--','-'};
hp = [];
for J = 1:2
    for q_i = -1:2:1
        if J == 1
            xi0 = [psi(z_i);0;q_i;hat_theta_i];
        else
            u0 = kappa_1(psi(z_i),[0;q_i;hat_theta_i]); %Fix initial input
            V2_0 = V_2(psi(z_i),[0;q_i;hat_theta_i;u0],theta);
            xi0 = [psi(z_i);0;q_i;hat_theta_i;u0;V2_0];
        end
        options = odeset('AbsTol',1e-8,'RelTol',1e-8);
        eval(['[t,j,xi] = HyEQsolver(@F_cl' num2str(J) ...
            ', @G_Omega' num2str(J) ...
            ', @C_Omega' num2str(J) ...
            ', @D_Omega' num2str(J) ...
            ', xi0, TSPAN, JSPAN, rule, options);']);
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
        plot(t,sqrt(sum((xi(:,6:7)-repmat(theta',[numel(t),1])).^2,2)),...
            'color',colors(J,:),'linestyle',linestyles{(q_i+3)/2})
        hold on
        figure(hf(3))
        if J == 1
            u = cell2mat(arrayfun(@(I) kappa_1(xi(I,1:3)',xi(I,4:7)'),1:numel(t),'UniformOutput',false));
            plot(t,sum(u.^2,1),'color',colors(J,:),'linestyle',linestyles{(q_i+3)/2})
        else
            u = xi(:,8:9)';
            plot(t,sum(u.^2,1),'color',colors(J,:),'linestyle',linestyles{(q_i+3)/2})
        end
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
text(z_0(1),z_0(2),'N','horizontalalignment','center','verticalalignment','middle')
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

figure(hf(3))
grid on
legend(hp,{'Section VI-B with q(0,0)=-1',...
    'Section VI-B with q(0,0)=1',...
    'Section VI-C with q(0,0)=-1',...
    'Section VI-C with q(0,0)=1'});
ylabel('|u(t)|')
xlabel('t')

% Problem Setup
function x = psi(z)
    global z_0 r
    x = [(z-z_0)/norm(z-z_0);log(norm(z-z_0)-r)];
end

function out = Dpsi(z)
    global z_0 r
    out= [PT((z-z_0)/norm(z-z_0))/norm(z-z_0);1/(norm(z-z_0)-r)*(z-z_0)'/norm(z-z_0)];
end

function out = DDpsi(z)
    global z_0 r
    s = (z-z_0)/norm(z-z_0);
    aux = Dpsi(z);
    out1 = -s*s'/(norm(z-z_0)-r)^2+PT(s)/(norm(z-z_0)-r)/norm(z-z_0);
    out2 = -vec(PT(s))*(z-z_0)'./norm(z-z_0)^3-(kron(s,eye(2))+kron(eye(2),s))*aux(1:2,:)/norm(z-z_0);
    out = [out2(1:2,:);out1(1,:);out2(3:4,:);out1(2,:)];
end

function z = psi_inv(x)
    global z_0 r
    z = z_0+x(1:2)*(exp(x(3))+r);
end

function out = Dpsi_inv(x)
    global r
    out = [[1;0]*(exp(x(3))+r),[0;1]*(exp(x(3))+r),x(1:2)*exp(x(3))];
end

function out = PT(x)
    n = numel(x);
    out = eye(n)-x*x';
end

function out = f(x,q)
    out = zeros(3,1);
end

function out = H(x,q)
    out = Dpsi(psi_inv(x));
end

function out = W(x,q)
    out = Dpsi(psi_inv(x));
end

function dx = F(x,q,u,theta)
    dx = f(x,q) + H(x,q)*u + W(x,q)*theta;
end

% The Nominally Synergistic Controller
function out = delta(x,q)
    out = 1;
end

function out = U(x,q)
    if q'*x(1:2) ~= 1
        out = 1;
    else
        out = 0;
    end
end

function out = phi(x,q)
    out = [(q(2)*x(1)-q(1)*x(2))/(1-q'*x(1:2)),x(3)]';
end

function out = Dqphi(x,q)
    d1 = -x(2)/(1-q'*x(1:2))+(q(2)*x(1)-q(1)*x(2))*x(1)/(1-q'*x(1:2))^2;
    d2 = x(1)/(1-q'*x(1:2))+(q(2)*x(1)-q(1)*x(2))*x(2)/(1-q'*x(1:2))^2;
    out = [d1 d2; 0 0];
end

function out = Dphi(x,q)
    out = zeros(2,3);
    out(2,3) = 1;
    out(1,1) = (q(2)-x(2))/(1-q'*x(1:2))^2;
    out(1,2) = -(q(1)-x(1))/(1-q'*x(1:2))^2;
end

function out = Dq_Dphi(x,q)
    out = zeros(6,2);
    out(1,1) = 2*(q(2)-x(2))*x(1)/(1-q'*x(1:2))^3;
    out(1,2) = (1-q'*x(1:2))^(-2)+2*(q(2)-x(2))*x(2)/(1-q'*x(1:2))^3;
    out(3,1) = -(1-q'*x(1:2))^(-2)-2*(q(1)-x(1))*x(1)/(1-q'*x(1:2))^3;
    out(3,2) = -2*(q(1)-x(1))*x(2)/(1-q'*x(1:2))^3;
end

function out = DDphi(x,q)
    out = zeros(6,3);
    out(1,1) = 2*q(1)*q(2)/(1-q'*x(1:2))^2+2*(q(2)*x(1)-q(1)*x(2))*q(1)^2/(1-q'*x(1:2))^3;
    out(1,2) = q(2)^2/(1-q'*x(1:2))^2-q(1)^2/(1-q'*x(1:2))^2+2*(q(2)*x(1)-q(1)*x(2))*q(1)*q(2)/(1-q'*x(1:2))^3;
    out(3,1) = out(1,2);
    out(3,2) = -2*q(1)*q(2)/(1-q'*x(1:2))^2+2*(q(2)*x(1)-q(1)*x(2))*q(2)^2/(1-q'*x(1:2))^3;
end

function out = V_0(x,q)
    if U(x,q)
        x0 = psi(zeros(2,1));
        out = 0.5*norm(phi(x,q)-phi(x0,q))^2;
    else
        out = Inf;
    end
end

function out =  DxV_0(x,q)
    out = (phi(x,q)-phi(psi(zeros(2,1)),q))'*Dphi(x,q);
end

function out = DqV_0(x,q)
    x0 = psi(zeros(2,1));
    out = (phi(x,q)-phi(x0,q))'*(Dqphi(x,q)-Dqphi(x0,q));
end

function out = Gamma(q)
    global gamma
    x0 = psi(zeros(2,1));
    out = max([0,gamma-x0(1:2)'*q]);
end

function dq = F_c(x,q)
    dq = -Gamma(q)*PT(q)*DqV_0(x,q)';
end

function out = kappa_0(x,q)
    global Gamma_0
    x0 = psi(zeros(2,1));
    out = -Gamma_0*Dpsi(psi_inv(x))'*Dphi(x,q)'*(phi(x,q)-phi(x0,q));
end

function out = Dqkappa_0(x,q)
    global Gamma_0
    out = -Gamma_0*Dpsi(psi_inv(x))'*( ...
        kron( (phi(x,q)-phi(psi(zeros(2,1)),q))' ,eye(3)) * K(2,3) * Dq_Dphi(x,q)...
        + Dphi(x,q)' * ( ...
            Dqphi(x,q) - Dqphi(psi(zeros(2,1)),q) ...
        ) ...
    );
end

function out = Dxkappa_0(x,q)
    global Gamma_0
    phi0 = phi(psi(zeros(2,1)),q);
    psi_u = Dpsi(psi_inv(x));
    DV = kron((phi(x,q)-phi0)',eye(3))*K(2,3)*DDphi(x,q)+Dphi(x,q)'*Dphi(x,q);
    Dpsi_u = DDpsi(psi_inv(x))*Dpsi_inv(x);
    out = -Gamma_0*(kron((phi(x,q)-phi0)'*Dphi(x,q),eye(2))*K(3,2)*Dpsi_u+psi_u'*DV);
end

function out = varrho_V0(x,q)
    global gamma
    x0 = psi(zeros(2,1));
    if x(1:2)'*x0(1:2) == -1
        %out = (randi(2)*2-3)*[x0(2);-x0(1)];
        out = [x0(2);-x0(1)]; % to ensure determinism
    else
        rho = -(x0(1:2)+x(1:2))/norm((x0(1:2)+x(1:2)));
        if rho'*x(2:3) > gamma
            out = gamma*x0(1:2)+sqrt(1-gamma^2)*PT(x(1:2))*rho;
        else
            out = rho;
        end
    end
end

function out = nu_V0(x,q)
    next_q = varrho_V0(x,q);
    out = V_0(x,next_q);
end

function out = mu_V0(x,q)
    out = V_0(x,q)-nu_V0(x,q);
end

function dxi = F_cl0(xi)
    x = xi(1:3);
    q = xi(4:5);
    theta = zeros(2,1); % Since this is the nominal case
    dx = F(x,q,kappa_0(x,q),theta);
    dq = F_c(x,q);
    dxi = [dx;dq];
end

function out = C_0(xi)
    x = xi(1:3);
    q = xi(4:5);
    if mu_V0(x,q) <= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function out = D_0(xi)
    x = xi(1:3);
    q = xi(4:5);
    if mu_V0(x,q) >= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function next_xi = G_cl0(xi)
    x = xi(1:3);
    q = xi(4:5);
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
    global hat_W
    q = x_c1(1:2);
    hat_theta = x_c1(3:4);
    u = kappa_0(x,q)-hat_W*hat_theta;
end

function out = Dxkappa_1(x,x_c1)
    q = x_c1(1:2);
    hat_theta = x_c1(3:4);
    out = Dxkappa_0(x,q)-kron(hat_theta',eye(2))*Dx_hat_W(x,q);
end

function out = Dqkappa_1(x,x_c1)
    q = x_c1(1:2);
    hat_theta = x_c1(3:4);
    out = Dqkappa_0(x,q)-kron(hat_theta',eye(2))*Dq_hat_W(x,q);
end

function out = Dx_hat_W(x,q)
    out = zeros(4,3); %hat_W = eye(2), so its derivative is 0
end
function out = Dq_hat_W(x,q)
    out = zeros(4,2); %hat_W = eye(2), so its derivative is 0
end

function out = V_1(x,x_c1,theta)
    global Gamma_1
    q = x_c1(1:2);
    hat_theta = x_c1(3:4);
    out = V_0(x,q)+0.5*norm(Gamma_1^(-1/2)*(theta-hat_theta))^2;
end

function dx_c1 = F_c1(x,x_c1)
    global Gamma_1 
    q = x_c1(1:2);
    hat_theta = x_c1(3:4);
    dq = F_c(x,q);
    dhat_theta = Gamma_1*Proj(Dpsi(psi_inv(x))'*DxV_0(x,q)',hat_theta);
    dx_c1 = [dq;dhat_theta];
end

function out = nu_V1(x,x_c1,theta)
    q = x_c1(1:2);
    out = nu_V0(x,q);
end

function out = varrho_V1(x,x_c1,theta)
    q = x_c1(1:2);
    hat_theta = x_c1(3:4);
    out = [varrho_V0(x,q);
        theta];
end

function out = mu_V1(x,x_c1,theta)
    global Gamma_1
    q = x_c1(1:2);
    hat_theta = x_c1(3:4);
    out = mu_V0(x,q)+0.5*norm(Gamma_1^{-1/2}*(theta-hat_theta))^2;
end

function out = min_mu_V1(x,x_c1)
    global theta_0
    q = x_c1(1:2);
    hat_theta = x_c1(3:4);
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
    q = x_c1(1:2);
    hat_theta = x_c1(3:4);
    out = [varrho_V0(x,q);
        hat_G(hat_theta)];
end

function out = C_Omega1(xi)
    x = xi(1:3);
    x_c1 = xi(4:7);
    q = x_c1(1:2);
    if min_mu_V1(x,x_c1) <= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function dxi = F_cl1(xi)
    global theta
    x = xi(1:3);
    x_c1 = xi(4:7);
    q = x_c1(1:2);
    dx = F(x,q,kappa_1(x,x_c1),theta);
    dx_c1 = F_c1(x,x_c1);
    dxi = [dx;dx_c1];
end

function out = D_Omega1(xi)
    x = xi(1:3);
    x_c1 = xi(4:7);
    q = x_c1(1:2);
    if min_mu_V1(x,x_c1) >= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function next_xi = G_Omega1(xi)
    x = xi(1:3);
    x_c1 = xi(4:7);
    next_x = x;
    next_x_c1 = G_c1(x,x_c1);
    next_xi = [next_x;
        next_x_c1];
end

% Backstepping
function out = v(x,x_c2)
    global Gamma_2
    x_c1 = x_c2(1:4);
    q = x_c1(1:2);
    u = x_c2(5:6);
    out = W(x,q)'*(DxV_0(x,q)'-Dxkappa_1(x,x_c1)'*Gamma_2^(-1)*(u-kappa_1(x,x_c1)));
end

function out = f_u(x,x_c2)
    global Gamma_1 Gamma_2 k_u hat_W
    x_c1 = x_c2(1:4);
    q = x_c1(1:2);
    hat_theta = x_c1(3:4);
    u = x_c2(5:6);
    out = -hat_W*Gamma_1*Proj(v(x,x_c2),hat_theta)-k_u*(u-kappa_1(x,x_c1))...
        -Gamma_2*H(x,q)'*DxV_0(x,q)'+Dxkappa_1(x,x_c1)*F(x,q,u,hat_theta);
end

function dx_c2 = F_c2(x,x_c2)
    global Gamma_1
    x_c1 = x_c2(1:4);
    q = x_c1(1:2);
    hat_theta = x_c1(3:4);
    u = x_c2(5:6);
    dx_c2 = [F_c(x,q);
        Gamma_1*Proj(v(x,x_c2),hat_theta);
        f_u(x,x_c2)+Dqkappa_1(x,x_c1)*F_c(x,q)];
end

function out = V_2(x,x_c2,theta)
    global Gamma_2
    x_c1 = x_c2(1:4);
    u = x_c2(5:6);
    out = V_1(x,x_c1,theta);
    if ~isinf(out)
        out = out+0.5*norm(Gamma_2^(-1/2)*(u-kappa_1(x,x_c1)))^2;
    end
end

function g_c2 = G_c2(x,x_c2)
    x_c1 = x_c2(1:4);
    g_c1 = G_c1(x,x_c1);
    g_u = kappa_1(x,g_c1);
    g_c2 = [g_c1;g_u];
end

function out = min_mu_V2(x,x_c2)
    global Gamma_2
    x_c1 = x_c2(1:4);
    u = x_c2(5:6);
    out = min_mu_V1(x,x_c1);
    if ~isinf(out)
        out = out+0.5*norm(Gamma_2^(-1/2)*(u-kappa_1(x,x_c1)))^2;
    end
end

function out = C_Omega2(xi)
    x = xi(1:3);
    x_c2 = xi(4:9);
    q = x_c2(1:2);
    if min_mu_V2(x,x_c2) <= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function dxi = F_cl2(xi)
    global theta Gamma_2 k_u
    x = xi(1:3);
    x_c2 = xi(4:9);
    x_c1 = x_c2(1:4);
    q = x_c2(1:2);
    hat_theta = x_c2(3:4);
    u = x_c2(5:6);
    dx = F(x,q,u,theta);
    dx_c2 = F_c2(x,x_c2);
    dV2 = DxV_0(x,q)*F(x,q,kappa_0(x,q),zeros(2,1))...
        -k_u*(u-kappa_1(x,x_c1))'*Gamma_2^(-1)*(u-kappa_1(x,x_c1));
    dxi = [dx;dx_c2;dV2];
end

function out = D_Omega2(xi)
    x = xi(1:3);
    x_c2 = xi(4:9);
    q = x_c2(1:2);
    if min_mu_V2(x,x_c2) >= delta(x,q)
        out = 1;
    else
        out = 0;
    end
end

function next_xi = G_Omega2(xi)
    global theta
    x = xi(1:3);
    x_c2 = xi(4:9);
    next_x = x;
    next_x_c2 = G_c2(x,x_c2);
    next_xi = [next_x;next_x_c2;V_2(next_x,next_x_c2,theta)];
end