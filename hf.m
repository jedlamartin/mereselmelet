clear; close all;

%% Becsléselméleti feladatok - 1. feladat

% Paraméterek
mu_A = 1;
mu_B = 2;
mu_C = 1;

sigma_a = 0.1;
rho = 0.2;

t0 = 10e-3;
f0 = 50;
N = 5;

fs = 250;
dt = 1/fs;
t = (t0 + (0:N-1)*dt)';

sigma_w = 0.2;

% Kovariancia mátrix
Caa = sigma_a^2 * [1 rho rho^2; rho 1 rho; rho^2 rho 1];
Cww = sigma_w^2 * eye(length(t));

% Paramétervektor generálása
mu_a = [mu_A; mu_B; mu_C];
a = mvnrnd(mu_a, Caa)';
A = a(1);
B = a(2);
C = a(3);

% Zajvektor generálása
w = sigma_w * randn(size(t));

% A jel generálása
z_clean = A*sin(2*pi*f0*t) + B*cos(2*pi*f0*t) + C;
z = z_clean + w;

%% 2. feladat

% Paraméterek felírása
U = [sin(2*pi*f0*t), cos(2*pi*f0*t), ones(length(t),1)];

% MS becslő
a_MS = mu_a + inv(U'*inv(Cww)*U + ...
    inv(Caa)) * U'*inv(Cww) * (z - U*mu_a);

% MS becslő torzítása
b_a_MS = (eye(3) - inv(U'*inv(Cww)*U + ...
    inv(Caa)) * U'*inv(Cww)*U) * mu_a + ...
    (inv(U'*inv(Cww)*U + inv(Caa)) * U'*inv(Cww)*U - eye(3)) * a;

% MS becslő kovarianciája
Caaz_MS = inv(U'*inv(Cww)*U + inv(Caa));

%% 3. feladat

% ML becslő
a_ML = inv(U'*inv(Cww)*U) * U'*inv(Cww)*z;

% ML becslő torzítása
b_a_ML = zeros(3,1);

% ML becslő kovarianciája
Caaz_ML = inv(U'*inv(Cww)*U);

%% 4. feladat

% LS becslő
a_LS = inv(U'*U) * U'*z;

% LS becslő torzítása
b_a_LS = zeros(3,1);

% LS becslő négyzetes hibája
Jaaz_LS = z'*(z-U*a_LS);

%% 5. feladat
% MS becslőből
D_hat_MS = sqrt(a_MS(1)^2 + a_MS(2)^2);
phi_hat_MS = atan2(a_MS(2), a_MS(1));

% Parciális deriváltak
dD_dA_MS = a_MS(1) / D_hat_MS;
dD_dB_MS = a_MS(2) / D_hat_MS;

dphi_dA_MS = -a_MS(2) / D_hat_MS^2;
dphi_dB_MS = a_MS(1) / D_hat_MS^2;

% Részmátrix Caa a becsült kovariancia mátrixból
C_AB_MS = Caaz_MS(1:2, 1:2);

var_D_MS = [dD_dA_MS, dD_dB_MS] * C_AB_MS * [dD_dA_MS; dD_dB_MS];
var_phi_MS = [dphi_dA_MS, dphi_dB_MS] * C_AB_MS * [dphi_dA_MS; dphi_dB_MS];

% ML becslőből
D_hat_ML = sqrt(a_ML(1)^2 + a_ML(2)^2);
phi_hat_ML = atan2(a_ML(2), a_ML(1));

% Parciális deriváltak
dD_dA_ML = a_ML(1) / D_hat_ML;
dD_dB_ML = a_ML(2) / D_hat_ML;

dphi_dA_ML = -a_ML(2) / D_hat_ML^2;
dphi_dB_ML = a_ML(1) / D_hat_ML^2;

% Részmátrix Caa a becsült kovariancia mátrixból
C_AB_ML = Caaz_ML(1:2, 1:2);

var_D_ML = [dD_dA_ML, dD_dB_ML] * C_AB_ML * [dD_dA_ML; dD_dB_ML];
var_phi_ML = [dphi_dA_ML, dphi_dB_ML] * C_AB_ML * [dphi_dA_ML; dphi_dB_ML];

% LS becslőből
D_hat_LS = sqrt(a_LS(1)^2 + a_LS(2)^2);
phi_hat_LS = atan2(a_LS(2), a_LS(1));

% Parciális deriváltak
dD_dA_LS = a_LS(1) / D_hat_LS;
dD_dB_LS = a_LS(2) / D_hat_LS;

dphi_dA_LS = -a_LS(2) / D_hat_LS^2;
dphi_dB_LS = a_LS(1) / D_hat_LS^2;

% % Részmátrix Caa a becsült kovariancia mátrixból
% C_AB_LS = Jaaz_LS(1:2, 1:2);
% 
% var_D_LS = [dD_dA_LS, dD_dB_LS] * C_AB_LS * [dD_dA_LS; dD_dB_LS];
% var_phi_LS = [dphi_dA_LS, dphi_dB_LS] * C_AB_LS * [dphi_dA_LS; dphi_dB_LS];


%% Multiszinuszos méréstechnika - 6. feladat

% Paraméterek megadása
M = 100;
bazis_meret=2*M+1;
N = 401;
n = 0:N-1;
x = ones(bazis_meret,1);

% Véletlen fázis, 0 a várható értéke
phi = 2*pi*rand(bazis_meret,1);

% Frekvenciák
theta = 2*pi*(0:bazis_meret-1)'/(bazis_meret);

% A jel generálása
c = exp(1i*(theta*n + repmat(phi, [1,N])));
u = sum(x.*c, 1);

% 0 fázis esetén
c_test = exp(1i*(theta*n));
u_test = sum(x.*c_test);

% Csúcsérték
u_max = max(abs(u));
u_test_max = max(abs(u_test));

%% 7. feladat

% g generálása
g = 1/bazis_meret * exp(-1i*(theta*n + repmat(phi, [1,N])));

% Az integrátorok kezdeti állapota
x_hat = zeros(bazis_meret,1);

% A kimenet
y_hat = zeros(1, length(u));

for n = 1:N
    y_hat(n) = sum(x_hat.*c(:,n));
    input = u(n) - y_hat(n);
    x_hat = x_hat + input*g(:,n);
end

figure;
plot(abs(y_hat-u), Color='black');
grid on;
xlabel('n');
ylabel('$$|\hat{y}-y|$$', 'Interpreter', 'latex');
fontsize(14,"points");


%% 8. feladat

% Paraméterek
r = 0.83;
system = @ D;

% Kimenet generálása
y = system(u, r);

% A kimenet
y_hat = zeros(1, length(u));

% Az integrátorok kezdeti állapota
x_hat = zeros(bazis_meret,1);

for n = 1:N
    y_hat(n) = sum(x_hat.*c(:,n));
    input = y(n) - y_hat(n);
    x_hat = x_hat + input*g(:,n);
end

th = 2*pi/length(x_hat)*(0:length(x_hat)-1);
figure;
plot(th/pi, abs(x_hat), Color='black');
grid on;
xlabel('$$\theta [\pi]$$', 'Interpreter','latex');
ylabel('$$|\hat{x}|$$', 'Interpreter', 'latex');
fontsize(14,"points");

figure;
plot(th/pi, angle(x_hat)/pi, Color='black');
grid on;
xlabel('$$\theta [\pi]$$', 'Interpreter','latex');
ylabel('$$angle(\hat{x}) [\pi]$$', 'Interpreter', 'latex');
fontsize(14,"points");

%% Modellillesztés, adaptív eljárások - 9. feladat

% Paraméterek
mu = 0.0001;
P = 15;

% Regressziós vektor
X = zeros(P,1);

% Súlytényezők
W = zeros(P,1);
W_buf = zeros(P, ceil(N/2));

y_reg = zeros(1, ceil(N/2));

for n = 1:ceil(N/2)
    y_reg(n) = W'*X;
    e = conj(y(n) - y_reg(n));
    W = W + 2*mu*e*X;
    X = [u(n); X(1:end-1)];
    W_buf(:,n) = W;
end


[~, W_maxI] = maxk(W, 5);
W_buf_top5 = abs(W_buf(W_maxI, :));

line_styles = {'-', '--', ':', '-.', '-'};
markers = {'o', 's', '^', 'd', 'x'};

figure;
hold on;
for i = 1:5
    plot(W_buf_top5(i,:), ...
        'LineStyle', line_styles{i}, ...
        'Marker', markers{i}, ...
        'MarkerIndices', 1:20:N/2, ...
        'Color', 'black');
end
hold off;
fontsize(14,"points");
xlabel('n');
ylabel('Amplitúdó');
legend("W("+string(W_maxI)+")", 'Location', 'best');
xlim([0, bazis_meret-1]);
grid on;


% r felére csökkentése
r = r/2;
y = system(u,r);

W_buf = zeros(P, floor(N/2));

y_reg = zeros(1, floor(N/2));

for n = 1:floor(N/2)
    y_reg(n) = W'*X;
    e = conj(y(n) - y_reg(n));
    W = W + 2*mu*e*X;
    X = [u(n+ceil(N/2)); X(1:end-1)];
    W_buf(:,n) = W;
end

W_buf_top5 = abs(W_buf(W_maxI, :));

line_styles = {'-', '--', ':', '-.', '-'};
markers = {'o', 's', '^', 'd', 'x'};

figure;
hold on;
for i = 1:5
    plot(W_buf_top5(i,:), ...
        'LineStyle', line_styles{i}, ...
        'Marker', markers{i}, ...
        'MarkerIndices', 1:20:N/2, ...
        'Color', 'black');
end
hold off;
fontsize(14,"points");
xlabel('n');
ylabel('Amplitúdó');
legend("W("+string(W_maxI)+")", 'Location', 'best');
xlim([0, bazis_meret-1]);
grid on;

%% 10. feladat
r = 0.83;
y = system(u,r);

mu = 0.001;
P = 7;

% Regressziós vektor
X = zeros(P,1);

% Súlytényezők
W = zeros(P,1);
W_buf = zeros(P, ceil(N/2));

y_reg = zeros(1, ceil(N/2));

for n = 1:ceil(N/2)
    y_reg(n) = W'*X;
    e = conj(y(n) - y_reg(n));
    W = W + 2*mu*e*X;
    X = [u(n); X(1:2); y(n); X(4:end-1)];
    W_buf(:,n) = W;
end

[~, W_maxI] = maxk(W, 5);
W_buf_top5 = abs(W_buf(W_maxI, :));


line_styles = {'-', '--', ':', '-.', '-'};
markers = {'o', 's', '^', 'd', 'x'};

figure;
hold on;
for i = 1:5
    plot(W_buf_top5(i,:), ...
        'LineStyle', line_styles{i}, ...
        'Marker', markers{i}, ...
        'MarkerIndices', 1:20:N/2, ...
        'Color', 'black');
end
hold off;
fontsize(14,"points");
xlabel('n');
ylabel('Amplitúdó');
legend("W("+string(W_maxI)+")", 'Location', 'best');
xlim([0, bazis_meret-1]);
grid on;

% r felére csökkentése
r = r/2;
y = system(u,r);

y_reg = zeros(1, floor(N/2));


for n = 1:floor(N/2)
    y_reg(n) = W'*X;
    e = conj(y(n) - y_reg(n));
    W = W + 2*mu*e*X;
    X = [u(n+ceil(N/2)); X(1:2); y(n+ceil(N/2)); X(4:end-1)];
    W_buf(:,n) = W;
end

W_buf_top5 = abs(W_buf(W_maxI, :));

line_styles = {'-', '--', ':', '-.', '-'};
markers = {'o', 's', '^', 'd', 'x'};

figure;
hold on;
for i = 1:5
    plot(W_buf_top5(i,:), ...
        'LineStyle', line_styles{i}, ...
        'Marker', markers{i}, ...
        'MarkerIndices', 1:20:N/2, ...
        'Color', 'black');
end
hold off;
fontsize(14,"points");
xlabel('n');
ylabel('Amplitúdó');
legend("W("+string(W_maxI)+")", 'Location', 'best');
xlim([0, bazis_meret-1]);
grid on;

%% 11. feladat
% Eredeti értékek visszaállítása
r = r*2;
y = system(u,r);

N_new = round(N/5);

sigma_n = u_max*0.02;
% sigma_w = 0;
sigma_w = 0.5;

% Kovariancia mártix
R = sigma_n^2;
Q = sigma_w^2 * eye(4);

% Az állapotváltozós leírás mátrixainak meghatározása
A = [0, 0, 0, r; 1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0];
B = [(1-r)/2; 0; 0; 0];
C = [0, 0, 1, 0];
D = 0;

x_hat = zeros(4,1);
P = eye(4);
trace_P = zeros(4,N);

for n = 1:ceil(N/2)
    e = y(n) - C * x_hat;

    % 2. Kalman-nyereség (G) frissítése
    G = A * P * C' * inv(C * P * C' + R);

    % 3. Állapotbecslés frissítése
    x_hat = A * x_hat + G * e;

    % 4. Hibakovariancia (P) frissítése
    P = (A - G * C) * P * A' + Q;

    trace_P(n) = trace(P);
end

figure;
plot(trace_P(1:N_new), Color='black');
grid on;
xlabel('Time step (n)');
ylabel('trace(P(n))');
fontsize(14,"points");

% r felére csökkentés
r = r/2;
y = system(u,r);

for n = 1:floor(N/2)
    e = y(n+ceil(N/2)) - C * x_hat;

    % 2. Kalman-nyereség (G) frissítése
    G = A * P * C' * inv(C * P * C' + R);

    % 3. Állapotbecslés frissítése
    x_hat = A * x_hat + G * e;

    % 4. Hibakovariancia (P) frissítése
    P = (A - G * C) * P * A' + Q;

    trace_P(n+ceil(N/2)) = trace(P);
end

figure;
plot(trace_P(ceil(N/2) : N/2+N_new), Color='black');
grid on;
xlabel('Time step (n)');
ylabel('trace(P(n))');
fontsize(14,"points");
ax = gca;
ax.YAxis.TickLabelFormat = '%.4f';