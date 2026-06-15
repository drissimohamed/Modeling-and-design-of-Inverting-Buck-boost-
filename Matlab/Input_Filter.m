%% ═══════════════════════════════════════════════════════════
%  Input Filter Analysis — Updated Design
%  Inverting Buck-Boost PCMC
%  Engineer: Mohamed Drissi
%% ═══════════════════════════════════════════════════════════

clear; clc; close all;

%% ── Converter Operating Point ────────────────────────────────
Vdc  = 24;
Vout = 12;
Iout = 3;
Pout = Vout * Iout;
Iin  = Pout / Vdc;
fs   = 300e3;
D    = Vout/(Vdc+Vout);
D1   = 1-D;
R    = Vout/Iout;
L    = 22e-6;
C    = 78.8e-6;
Ri   = 0.3;
rc   = 0.002;
rl   = 50e-3;
Vref = 1.215;
Tsw  = 1/fs;

fprintf('Operating Point\n');
fprintf('  D  = %.4f\n', D);
fprintf('  D1 = %.4f\n', D1);
fprintf('  R  = %.2f Ohm\n', R);

%% ── Input Filter Parameters (NEW VALUES) ─────────────────────
Lf       = 10e-6;     % Filter inductor [H] — sized for fc=2kHz with Cf_derated
DCR_Lf   = 18e-3;      % Inductor DCR [Ohm]
Cf       = 5.45e-6;    % Cin derated @ 24V = 
ESR_Cf   = 2e-3;       % Cin effective ESR (2 caps in //) [Ohm]
Rdamp    = 0.4;        % Damping resistor [Ohm]
Cdamp    = 220e-6;     % Damping capacitor — polymer electrolytic [F]
ESR_Cd   = 11e-3;      % Cdamp ESR (polymer electrolytic) [Ohm]
Cin      = 53.2e-6;    % Input Capacitance [F]
ESR_Cin  = 2e-3;
fprintf('\nFilter Parameters\n');
fprintf('  Lf       = %.0f uH\n',  Lf*1e6);
fprintf('  DCR_Lf   = %.0f mOhm\n', DCR_Lf*1e3);
fprintf('  Cf       = %.1f uF (derated @ 24V)\n', Cf*1e6);
fprintf('  ESR_Cf   = %.0f mOhm\n', ESR_Cf*1e3);
fprintf('  Rdamp    = %.1f Ohm\n',  Rdamp);
fprintf('  Cdamp    = %.0f uF\n',   Cdamp*1e6);
fprintf('  ESR_Cd   = %.0f mOhm\n', ESR_Cd*1e3);

%% ── Key Filter Frequencies ───────────────────────────────────
fc_filter = 1/(2*pi*sqrt(Lf*Cf));
Z0        = sqrt(Lf/Cf);
Q_filter  = Rdamp * sqrt(Cf/Lf);
A_dB      = 40*log10(fs/fc_filter);

fprintf('\nFilter Key Values\n');
fprintf('  fc_filter = %.2f kHz\n', fc_filter/1e3);
fprintf('  Z0        = %.3f Ohm\n', Z0);
fprintf('  Q         = %.3f\n',     Q_filter);
fprintf('  A @ fs    = %.1f dB\n',  A_dB);

%% ── TF Objects ───────────────────────────────────────────────
s = tf('s');

%% ── Filter Output Impedance WITH PARASITICS ──────────────────

% Inductor: Lf + DCR in series
Z_Lf = s*Lf + DCR_Lf;

% Cf with ESR
Z_Cf = 1/(s*Cf) + ESR_Cf;

%Cin with ESR
Z_Cin = 1/(s*Cin) + ESR_Cin;

% Damping branch: Rdamp + Cdamp(with ESR) in series
Z_damp = Rdamp + ESR_Cd + 1/(s*Cdamp);

% Cin parallel with damping branch
Z_Cin_parallel = (Z_Cin * Z_damp) / (Z_Cin + Z_damp);

% L in series with Cf
Z_Cf_series = Z_Lf ;

% Filter output impedance (Lf in series with Cf in //  with parallel combination)
Z_filter_total = (Z_Cf_series * Z_Cin_parallel) / (Z_Cf_series + Z_Cin_parallel);

%% ── Converter Input Impedance ────────────────────────────────
Sn  = (Vdc/L)*Ri;
Sf  = (Vout/L)*Ri;
Se  = 90000;
mc  = 1 + (Se/Sn);

Vap = Vdc + Vout;
Vcp = Vout;
Vc  = 1.7;
Ic  = (Vc/Ri) - D*Tsw*Se - Vcp*(1-D)*Tsw/(2*L);

D   = Vcp/Vap;
D1  = 1-D;
ki  = D/Ri;
go  = (Tsw/L)*(D1*Se/Sn + 0.5 - D);
gf  = D*go - (D*D1)*Tsw/(2*L);
gi  = D*(gf - Ic/Vap);
gr  = (Ic/Vap) - go*D;
ko  = 1/Ri;
C3  = 1/(L*(fs*pi)^2);
req = 1/(go + 1/R + 1/rc);
Ro  = ((gi-gf+go+gr)*R+1) / ((gf*gr+gi*go)*R+gi);

Tau = (L/R)*fs;

tau1  = L / (((R*(gf*gr+gi*go)+gi)/(gf*gr+gi*go)) + rl);
tau2  = (((R*gi)/(R*(gf*gr+gi*go)+gi)) + rc) * C;
tau3  = ((R*gi)/(R*(gf*gr+gi*go)+gi)) * C3;
b1    = tau1 + tau2 + tau3;

tau12 = (rc+R)*C;
tau13 = 1/((gf*gr/gi) + go) * C3;
tau23 = ((req*gi)/(gi*(2*req*gf+1)+req*gf*gr)) * C3;
b2    = tau1*tau12 + tau1*tau13 + tau2*tau23;
b3    = tau1*tau12*tau13;

wz1 = (((1-D)^3/(2*Tau))*(1+2*Se/Sn)+1+D) / (R*C);
Qp  = 1/(pi*(mc*D1-0.5));
wn  = pi/Tsw;

Zin = Ro * tf([(1/wz1) 1],[b3 b2 b1 1]) ...
         * tf([(1/wn^2) (1/(wn*Qp)) 1], 1);

fprintf('\nConverter Input Impedance\n');
fprintf('  Ro = %.4f Ohm\n', Ro);
fprintf('  Qp = %.4f\n',     Qp);

%% Plotting
w = logspace(2,6,1000);
figure;
bodemag(Z_filter_total,'b',Zin,'r',w);
grid on;
legend('Filter Output Impedance', ...
       'IBB Input Impedance', ...
       'Location','southwest');
title('Input Filter Output Impedance vs IBB Input Impedance')