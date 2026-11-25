```matlab
close all
clear all
clc

%% Define parameters
ky = 5 * 10^7; % (N/m)
alpha = 0.1;
beta = 2 * 10^(-4);

%% Import the structure and assemble the matrices
[file_i,xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr] = loadstructure;

[M,K] = assem(incidenze,l,m,EA,EJ,gamma,idb);
C = alpha*M + beta*K;

% Lumped spring
i_dof55 = idb(55, :);
K(i_dof55(2), i_dof55(2)) = K(i_dof55(2), i_dof55(2)) + ky;

dis_stru(posiz,l,gamma,xy,pr,idb,ndof);

% Define the constrained nodes
nConstr = [1 10 18 23];
num_nConstr = length(nConstr);
ndof_nConstr = num_nConstr*3;

% Matrices associated with free degrees of freedom
Mff = M(1:end-ndof_nConstr, 1:end-ndof_nConstr);
Kff = K(1:end-ndof_nConstr, 1:end-ndof_nConstr);
Cff = C(1:end-ndof_nConstr, 1:end-ndof_nConstr);

%% Eigenvalue analysis
[Phi,D] = eig(Mff \ Kff);

omega2 = diag(D);
omega = sqrt(omega2);
freq = omega/(2*pi);

% Sort natural frequencies
[freq_sorted,ind] = sort(freq);
Phi = Phi(:,ind);

% Select modes between 0 and 6 Hz
selected = find(freq_sorted > 0 & freq_sorted <= 6);
selected = selected(1:4);

freq_selected = freq_sorted(selected);

% Display selected natural frequencies
disp('Eigenfrequencies between 0 and 6 Hz:')
disp(freq_selected)

phi_4modes = zeros(length(Phi), length(selected));
scale_factor = 10;

for k = 1:length(selected)

    mode = Phi(:,selected(k));
    phi_4modes(:, k) = mode;

    figure;
    title(['Mode ', num2str(k)], 'Interpreter', 'latex', 'FontSize', 14);
    subtitle(['$f$ = ', num2str(freq_selected(k)), ' Hz'], 'Interpreter', 'latex', 'FontSize', 12);
    diseg2(mode,scale_factor,incidenze,l,gamma,posiz,idb,xy);
end

%% Point 3: Horizontal force applied at A, output at A and B

% External force
F0 = zeros(ndof,1);
F0(idb(9, 1)) = 1;
om = (0:0.01:6)*2*pi;

for ii = 1:length(om)
    A = -om(ii)^2*Mff + 1i*om(ii)*Cff + Kff;
    X(:,ii) = A\F0;
end

% FRF: input at A, horizontal output at A
figure;
sgtitle(['FRF: input in A, horizontal output in A'], "Interpreter", "latex", 'FontSize', 14, 'FontWeight', 'bold');

subplot(2, 1, 1);
semilogy(om/(2*pi), abs(X(idb(9, 1), :)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Module (m/N)', 'Interpreter','latex');
hold on
xline(freq_selected, '--');
grid on

subplot(2, 1, 2);
plot(om/(2*pi), angle(X(idb(9, 1), :)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Phase (rad)', 'Interpreter','latex');
hold on
xline(freq_selected, '--');
grid on

% FRF: input at A, vertical output at B
figure;
sgtitle(['FRF: input in A, vertical output in B'], "Interpreter", "latex", 'FontSize', 14);

subplot(2, 1, 1);
semilogy(om/(2*pi), abs(X(idb(55, 2), :)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Module (m/N)', 'Interpreter','latex');
hold on
xline(freq_selected, '--');
grid on

subplot(2, 1, 2);
plot(om/(2*pi), angle(X(idb(55, 2), :)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Phase (rad)', 'Interpreter','latex');
hold on
xline(freq_selected, '--');
grid on

%% Point 4: Vertical force at B, horizontal output at A
% FEM frequency response
F0 = zeros(ndof,1);
F0(idb(55, 2)) = 1;
om = (0:0.01:6)*2*pi;

for ii = 1:length(om)
    A = -om(ii)^2*Mff + 1i(om(ii)*Cff) + Kff;
    X(:,ii) = A\F0;
end

figure;
sgtitle(['FRF: input in B, horizontal output in A'], 'Interpreter', 'latex', 'FontSize', 14);

subplot(2, 1, 1);
h1 = semilogy(om/(2*pi), abs(X(idb(9, 1), :)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Module (m/N)', 'Interpreter','latex');
grid on;
xline(freq_selected, '--');
hold on;

subplot(2, 1, 2);
plot(om/(2*pi), angle(X(idb(9, 1), :)), 'LineWidth',1.5);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Phase (rad)', 'Interpreter','latex');
grid on
hold on;

% Modal matrices
% Modal superposition using the first two mode shapes
ii = 1:2;
Phi = phi_4modes(:,ii);

Mmod = Phi'*Mff*Phi;
Kmod = Phi'*Kff*Phi;
Cmod = Phi'*Cff*Phi;
Fmod = Phi'*F0;

for ii = 1:length(om)
    xx_mod(:,ii) = (-om(ii)^2*Mmod+1i(om(ii)*Cmod)+Kmod)\Fmod;
end

xx_m = Phi * xx_mod;

subplot(2, 1, 1);
h2 = semilogy(om/(2*pi), abs(xx_m(idb(9, 1), :)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Module (m/N)', 'Interpreter','latex');
xline(freq_selected, '--');
grid on;
legend([h1, h2], 'FEM', 'Modal');

subplot(2, 1, 2);
plot(om/(2*pi), angle(xx_m(idb(9, 1), :)), 'LineWidth',1.5);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Phase (rad)', 'Interpreter','latex');
xline(freq_selected, '--');
grid on
legend('FEM', 'Modal');
hold off;

%% Point 5: Static response to snow loading

p0 = -7000; % [N/m]
p0G = [0 p0]';

FG = zeros(3*nnod,1);

for ii = 55:65
    gammaii = gamma(ii);
    Lk = l(ii);

    p0L = [cos(gammaii) sin(gammaii); -sin(gammaii) cos(gammaii)]*p0G;

    FkL = [Lk/2; 0; 0; Lk/2; 0; 0]*p0L(1) + ...
          [0; Lk/2; Lk^2/12; 0; Lk/2; -Lk^2/12]*p0L(2);

    FkG = [cos(gammaii) -sin(gammaii) 0 0 0 0;
           sin(gammaii) cos(gammaii) 0 0 0 0;
           0 0 1 0 0 0;
           0 0 0 cos(gammaii) -sin(gammaii) 0;
           0 0 0 sin(gammaii) cos(gammaii) 0;
           0 0 0 0 0 1]*FkL;

    Ek = zeros(6,3*nnod);

    for jj = 1:6
        Ek(jj,incidenze(ii,jj)) = 1;
    end

    FG = FG + Ek'*FkG;
end

F0 = FG(1:ndof);

for ii = 1:length(om)
    A = -om(ii)^2*Mff + 1i(om(ii)*Cff) + Kff;
    X(:,ii) = A\F0;
end

idof_YA = idb(4,2);
FRF_YA = X(idof_YA,:);

nodi_peso = [9 48 17 49 50 22 51 25 52 53 54 55];

xF = Kff\F0;

figure;
title(['Static deflection due to snow loading'], 'Interpreter', 'latex', 'FontSize', 14);
xlabel(['x [m]'],'Interpreter','latex');
ylabel(['y [m]'], 'Interpreter','latex');
diseg2(xF,10,incidenze,l,gamma,posiz,idb,xy);

% Maximum displacement
diff = zeros(nnod, 1);

for k = 1:nnod
    if idb(k, 2) <= ndof
        diff(k) = abs(xF(idb(k, 2)));
    end
end

[max_disp, loc] = max(diff);

%% Point 6: Moving load analysis

g = 9.81;       % Gravitational acceleration [m/s^2]
m = 90;         % Moving mass [kg]
v0 = 2;         % Initial velocity [m/s]
a1 = 3.5;       % First acceleration phase [m/s^2]
a2 = 1.5;       % Second acceleration phase [m/s^2]
vf1 = 23.728;   % Reference velocity [m/s]
L_arc = 63.65;  % Arc length [m]
L_ret = 79.85;  % Return length [m]

P = m*g;        % Moving load [N]

dt = 0.001;
time = 0:dt:30;
T1 = 6.208;
T2 = T1 + 2.487;

% Calculate velocity and displacement as functions of time
for i = 1:length(time)

    if time(i) <= T1

        if i == 1
            v(i) = v0 + a1*dt;
            s(i) = v(i)*dt + 1/2 * a1 * dt^2;
        else
            v(i) = v(i-1) + a1*dt;
            s(i) = s(i-1) + v(i-1)*dt + 1/2 * a1 * dt^2;
        end

    else
        v(i) = v(i-1) + a2*dt;
        s(i) = s(i-1) + v(i-1)*dt + 1/2 * a2 * dt^2;
    end
end

% Plot velocity and displacement
figure;

subplot(2,1,1)
plot(time, v);
hold on;
xline(T1, '--');
hold on;
xline(T2, '--');
grid on;

subplot(2,1,2)
plot(time, s);
hold on;
xline(T1, '--');
hold on;
xline(T2, '--');
grid on;

%% Map displacement s(t) onto the FE mesh

rampa = 55:65;
L_rampa = l(rampa);
L_prog = cumsum(L_rampa);

position_force = zeros(1, length(time));
local_xi = zeros(1, length(time));

for i = 1:length(time)

    % Identify the element containing the moving load
    idx = find(s(i) <= L_prog, 1, 'first');

    if ~isempty(idx)

        sel_el = rampa(idx);
        position_force(i) = sel_el;

        % Calculate the local coordinate within the element
        if idx == 1
            local_xi(i) = s(i);
        else
            local_xi(i) = s(i) - L_prog(idx-1);
        end

        % Calculate the distance from the second node
        local_b(i) = l(sel_el) - local_xi(i);

    else

        % The moving load has passed the end of the track
        position_force(i) = 0;

    end
end

local_a = local_xi(1:length(local_b));

for i = 1:length(local_a)

    f_y = -P*cos(gamma(position_force(i)));
    f_x = P*sin(gamma(position_force(i)));

    f1y(i) = (f_y*local_a(i)*local_b(i)*(local_b(i) - local_a(i))) / l(position_force(i))^3 + f_y*local_b(i)/l(position_force(i));
    f2y(i) = (f_y*local_a(i)*local_b(i)*(local_a(i) - local_b(i))) / l(position_force(i))^3 + f_y*local_a(i)/l(position_force(i));
    f1x(i) = f_x*(l(position_force(i)) - local_a(i))/l(position_force(i));
    f2x(i) = f_x*local_a(i)/l(position_force(i));

    m1(i) = f_y*local_a(i)*(local_b(i)^2) / (l(position_force(i))^2);
    m2(i) = -f_y*local_b(i)*(local_a(i)^2) / (l(position_force(i))^2);
end

% Transform force components from local to global coordinates
for ii = 1:length(local_a)

    gammaii = gamma(position_force(ii));

    f1_temp = [f1x(ii) f1y(ii)]';
    f2_temp = [f2x(ii) f2y(ii)]';

    f1_glob(:, ii) = [cos(gammaii) -sin(gammaii);
                      sin(gammaii) cos(gammaii)]*f1_temp;

    f2_glob(:, ii) = [cos(gammaii) -sin(gammaii);
                      sin(gammaii) cos(gammaii)]*f2_temp;
end

rogoredo = [f1_glob;
            m1;
            f2_glob;
            m2];

rogoredo_full = zeros(ndof, length(local_a));

for i = 1:length(local_a)
    posizione = incidenze(position_force(i), :);
    rogoredo_full(posizione, i) = rogoredo(:, i);
end

% Define the generalized forces for the modal equations
Q = phi_4modes' * rogoredo_full;
Q = [Q zeros(size(phi_4modes, 2), length(time) - length(local_a))];

M_modal = phi_4modes' * Mff * phi_4modes;
K_modal = phi_4modes' * Kff * phi_4modes;
C_modal = phi_4modes' * Cff * phi_4modes;

%% Generalized force plots

figure;
plot(time(1:length(local_a)), rogoredo_full(idb(50, 1), :), 'b-o');
hold on;
plot(time(1:length(local_a)), rogoredo_full(idb(50, 2), :), 'r-o');
hold on;
plot(time(1:length(local_a)), rogoredo_full(idb(50, 3), :), 'g-o');

figure;
plot(time(1:length(local_a)), rogoredo_full(incidenze(55:65, 2), :));

%% Modal equations of motion

% Define the time-dependent generalized force
Q_interp = @(t) interp1(time', Q', t, 'linear', 'extrap')';

% Number of retained modes
n = size(M_modal,1);

% Equations of motion
odefun = @(t, z) [z(n+1:end); M_modal \ (Q_interp(t) - C_modal*z(n+1:end) - K_modal*z(1:n))];

% Initial conditions
q0 = zeros(n, 1);
q_punto0 = zeros(n, 1);
z0 = [q0; q_punto0];

% Time interval
tspan = [time(1) time(end)];

% Numerical solution
[t, z] = ode45(odefun, tspan, z0);

% Recover modal coordinates and velocities
q = z(:,1:n);
q_punto = z(:,n+1:end);

% Plot modal coordinates
figure
plot(t,q,'LineWidth',1.5)
title(['Modal coordinates'], 'Interpreter','latex', 'FontSize', 14);
legend('q1', 'q2', 'q3', 'q4')
xlabel(['time [s]'], 'Interpreter','latex', 'FontSize', 14);
ylabel(['q'], 'Interpreter','latex', 'FontSize', 14);
grid on

x1 = xline(T2, '--');
x1.Annotation.LegendInformation.IconDisplayStyle = 'off';

x2 = xline(T1, '--');
x2.Annotation.LegendInformation.IconDisplayStyle = 'off';

%% Displacement of point B as a function of time

q = q';

% Recover physical displacements
x_phys = phi_4modes * q;

nodo_B = 55;
idx_By = idb(nodo_B, 2);
y_B_t = x_phys(idx_By, :);

%% Displacement response at point B

figure
plot(t, y_B_t, 'LineWidth', 1.5);
title(['Displacement of point B vs time'], 'Interpreter','latex', 'FontSize', 14);
xlabel(['time [s]'], 'Interpreter','latex', 'FontSize', 14)
ylabel(['Displacement [m]'], 'Interpreter','latex', 'FontSize', 14)
xline(T2, '--');
hold on
xline(T1, '--');
grid on
```