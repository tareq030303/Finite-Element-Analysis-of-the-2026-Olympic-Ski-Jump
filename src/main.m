close all
clear all
clc

%% Define parameters
ky = 5 * 10^7; % (N/m)
alpha = 0.1;
beta = 2 * 10^(-4);

%% Import the structure and assembly the matrices
[file_i,xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure;

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

% Matrices related to free nodes
Mff = M(1:end-ndof_nConstr, 1:end-ndof_nConstr);
Kff = K(1:end-ndof_nConstr, 1:end-ndof_nConstr);
Cff = C(1:end-ndof_nConstr, 1:end-ndof_nConstr);

%% Eigenvalue structure
[Phi,D] = eig(Mff \ Kff);

omega2 = diag(D);
omega  = sqrt(omega2);
freq = omega/(2*pi);

% SORT FREQUENCIES
[freq_sorted,ind] = sort(freq);
Phi = Phi(:,ind);

% KEEP MODES BETWEEN 0 AND 6 Hz
selected = find(freq_sorted > 0 & freq_sorted <= 6);
selected = selected(1:4);

freq_selected = freq_sorted(selected);

% DISPLAY RESULTS
disp('Eigenfrequencies between 0 and 6 Hz:')
disp(freq_selected)

phi_4modes = zeros(length(Phi), length(selected));
scale_factor = 10;

for k = 1:length(selected) 

    mode = Phi(:,selected(k));
    phi_4modes(:, k) = mode;     % note on dimension: 153 is the num of free dofs (55x3 - 12)

    figure;
    title(['Mode ', num2str(k)], 'Interpreter', 'latex', 'FontSize', 14);
    subtitle(['$f$ = ', num2str(freq_selected(k)), ' Hz'], 'Interpreter', 'latex', 'FontSize', 12);
    diseg2(mode,scale_factor,incidenze,l,gamma,posiz,idb,xy);
end

%% Punto 3: Horizontal force applied in A, output in A and B
% External force
F0 = zeros(ndof,1);
F0(idb(9, 1)) = 1;
om = (0:0.01:6)*2*pi;
for ii=1:length(om)
    A = -om(ii)^2*Mff + 1i*om(ii)*Cff + Kff;
    X(:,ii) = A\F0;
end

% plot input A output A
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

%plot input A output B
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


% % External force (negative) - it simply shows opposite phase
% F0 = zeros(ndof,1);
% F0(idb(9, 1)) = -1;
% om = (0:0.01:6)*2*pi;
% for ii=1:length(om)
%     A = -om(ii)^2*Mff + 1i*om(ii)*Cff + Kff;
%     X(:,ii) = A\F0;
% end
% 
% figure;
% subplot(2, 1, 1);
% semilogy(om/(2*pi), abs(X(idb(9, 1), :)));
% grid on
% subplot(2, 1, 2);
% plot(om/(2*pi), angle(X(idb(9, 1), :)));
% grid on 
% title('FRF fem NEGATIVE F: input A , horizontal output A');

%% Punto 4: Vertical force in B, horizontal output in A - comparison FEM and modal
% External force - FEM
F0 = zeros(ndof,1);
F0(idb(55, 2)) = 1;
om = (0:0.01:6)*2*pi;
for ii=1:length(om)
    A = -om(ii)^2*Mff + 1i*om(ii)*Cff + Kff;
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
% first 2 mode shapes, FRF in modal superposition approach
ii = 1:2;
Phi = phi_4modes(:,ii);
Mmod = Phi'*Mff*Phi;
Kmod = Phi'*Kff*Phi;
Cmod = Phi'*Cff*Phi;
Fmod = Phi'*F0;

for ii = 1:length(om)
    xx_mod(:,ii) = (-om(ii)^2*Mmod+1i*om(ii)*Cmod+Kmod)\Fmod;
end
xx_m = Phi * xx_mod;                   % notice: xx_mod is what we usually call q 

subplot(2, 1, 1);
h2 = semilogy(om/(2*pi), abs(xx_m(idb(9, 1), :)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Module (m/N)', 'Interpreter','latex');
xline(freq_selected, '--');
grid on;
legend([h1, h2], 'FEM', 'Modal');
subplot(2, 1, 2);
plot(om/(2*pi), angle(xx_m(idb(9, 1), :)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)', 'Interpreter','latex');
ylabel('Phase (rad)', 'Interpreter','latex');
xline(freq_selected, '--');
grid on 
legend('FEM', 'Modal');
hold off;


%% Punto 5
p0 = -7000; %[N/m]
p0G = [0 p0]';
FG = zeros(3*nnod,1);
for ii = 55:65 % cicle over elements involved
    gammaii = gamma(ii);
    Lk = l(ii);
    p0L = [cos(gammaii) sin(gammaii); -sin(gammaii) cos(gammaii)]*p0G;
    FkL = [Lk/2; 0 ; 0 ; Lk/2; 0 ; 0 ]*p0L(1) + [0 ; Lk/2; Lk^2/12; 0 ; Lk/2; -Lk^2/12]*p0L(2);
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
for ii=1:length(om)
    A = -om(ii)^2*Mff + 1i*om(ii)*Cff + Kff;
    X(:,ii)=A\F0;
end
idof_YA = idb(4,2);
FRF_YA = X(idof_YA,:);

nodi_peso = [9 48 17 49 50 22 51 25 52 53 54 55];

xF = Kff\F0;
figure;
title(['Static deflection due to the snow'], 'Interpreter', 'latex', 'FontSize', 14);
xlabel(['x [m]'],'Interpreter','latex');
ylabel(['y [m]'], 'Interpreter','latex');
diseg2(xF,10,incidenze,l,gamma,posiz,idb,xy);

% Massimo spostamento
diff = zeros(nnod, 1);
for k = 1:nnod
    if idb(k, 2) <= ndof
        diff(k) = abs(xF(idb(k, 2)));
    end
end
[max_disp, loc] = max(diff);


%% Punto 6
% Dati
g = 9.81;       % m/s^2
m = 90;         % kg
v0 = 2;         % m/s
a1 = 3.5;       % m/s^2
a2 = 1.5;       % m/s^2
vf1 = 23.728;   % m/s
L_arc = 63.65;  % m
L_ret = 79.85;  % m

P = m*g;        % N; moving load

dt = 0.001;
time = 0:dt:30;
T1 = 6.208;
T2 = T1 + 2.487;

% Calculation of the velocity and displacement in function of time
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


% v_plot = [v(1:length(local_a)) zeros(1, length(position_force) - length(local_a))];
% s_plot = [s(1:length(local_a)) zeros(1, length(position_force) - length(local_a))];
% Plots
figure;
subplot(2,1,1)
plot(time, v); %_plot);
hold on;
xline(T1, '--');
hold on;
xline(T2, '--');
grid on;
subplot(2,1,2)
plot(time, s); %_plot, "r");
hold on;
xline(T1, '--');
hold on;
xline(T2, '--');
grid on;

%% Mappatura s(t) sulla Mesh FE
rampa = 55:65;           % Indici degli elementi della pista
L_rampa = l(rampa);      % Vettore delle lunghezze degli elementi coinvolti
L_prog = cumsum(L_rampa); % Lunghezza progressiva (distanza di ogni nodo da A)

position_force = zeros(1, length(time)); % Pre-allocazione per velocità
local_xi = zeros(1, length(time));       % Coordinata locale a (distanza dal primo nodo dell'elemento)

for i = 1:length(time)
    % Trova il primo elemento k la cui fine è oltre la posizione attuale s(i)
    idx = find(s(i) <= L_prog, 1, 'first');
    
    if ~isempty(idx)
        sel_el = rampa(idx);          % Identifica l'elemento reale nella mesh globale
        position_force(i) = sel_el;
        
        % Calcolo della coordinata locale xi (distanza dal nodo iniziale dell'elemento k)
        if idx == 1
            local_xi(i) = s(i);
        else
            local_xi(i) = s(i) - L_prog(idx-1);
        end
        % Calcolo di b (distanza dal secondo nodo)
        local_b(i) = l(sel_el) - local_xi(i); 
    else

        % Il saltatore è fuori dalla pista (ha superato il punto B)
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

% Changing the reference system from local to global
for ii = 1:length(local_a)     % cycle over elements involved
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

% Define the components associated to the modal equation
Q = phi_4modes' * rogoredo_full;
Q = [Q zeros(size(phi_4modes, 2), length(time) - length(local_a))];
M_modal = phi_4modes' * Mff * phi_4modes;
K_modal = phi_4modes' * Kff * phi_4modes;
C_modal = phi_4modes' * Cff * phi_4modes;

%% Plot
figure;
plot(time(1:length(local_a)), rogoredo_full(idb(50, 1), :), 'b-o');
hold on;
plot(time(1:length(local_a)), rogoredo_full(idb(50, 2), :), 'r-o');
hold on;
plot(time(1:length(local_a)), rogoredo_full(idb(50, 3), :), 'g-o');


%%
figure;
plot(time(1:length(local_a)), rogoredo_full(incidenze(55:65, 2), :));


%% Mq'' + Cq' + Kq = Q(t)
% Define the time over which Q acts
% t_force = linspace(0, 10, length(local_a));

% Write Q in function of time
Q_interp = @(t) interp1(time', Q', t, 'linear', 'extrap')';

% Number of modes
n = size(M_modal,1);

% EoM
odefun = @(t, z) [z(n+1:end); M_modal \ (Q_interp(t) - C_modal*z(n+1:end) - K_modal*z(1:n))];

% Initial conditions
q0 = zeros(n, 1);
q_punto0 = zeros(n, 1);
z0 = [q0; q_punto0];

% Intervallo temporale
tspan = [time(1) time(end)];

% Soluzione
[t, z] = ode45(odefun, tspan, z0);

% Back to space and speed
q = z(:,1:n);
q_punto = z(:,n+1:end);

% Plot
figure
plot(t,q,'LineWidth',1.5)
title(['Modal coordinates'], 'Interpreter','latex', 'FontSize', 14);
%legend(q_plot)
legend('q1', 'q2', 'q3', 'q4')
xlabel(['time [s]'], 'Interpreter','latex', 'FontSize', 14);
ylabel(['q'], 'Interpreter','latex', 'FontSize', 14);
grid on
x1 = xline(T2, '--');
x1.Annotation.LegendInformation.IconDisplayStyle = 'off';
x2 = xline(T1, '--');
x2.Annotation.LegendInformation.IconDisplayStyle = 'off';
%% Displacement of point B vs time
q = q';

% Define x 
x_phys = phi_4modes * q;

nodo_B = 55;             % Sostituisci con il numero del tuo nodo B
idx_By = idb(nodo_B, 2); % indica il DOF verticale (y)
y_B_t = x_phys(idx_By, :);

% Qlag = phi_4modes * Q;

%% PLot

figure
plot(t, y_B_t, 'LineWidth', 1.5);
title(['Displacement of point B vs time'], 'Interpreter','latex', 'FontSize', 14);
xlabel(['time [s]'], 'Interpreter','latex', 'FontSize', 14)
ylabel(['space [m]'], 'Interpreter','latex', 'FontSize', 14)
xline(T2, '--');
hold on
xline(T1, '--');
grid on