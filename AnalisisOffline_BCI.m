% DETECCIÓN DE ONDA MU realtime(ERD/ERS) - ADS1299
% EMG - 1 Canal diferencial (M1 - M2)
% EEG - (C3 - [C3'+C3'']/2) y (C4 - [C4'+C4'']/2)
% 0.5 Segundos de calibración para disminución de respuesta transitoria en filtrado
% ventanas de referencia base de 5 segundos
% tiempo muerto entre referencia y actual (gap = 1 segundo)
% ventana de análisis actual 1 segundo
% desplazamiento cada 0.25 segundo
% persistencia ERD de 2 desplazamientos
% detección de picos gigantes como artefacto
% motor se mueve 1 segundo o ERS detiene antes
% No se calcula promedio si hubo artefactos
% =========================================================================

% ---------------------------------------------------------
% 1. Configuración inicial y lectura del archivo
% ---------------------------------------------------------

clear; close all; clc

filename = 'TEST_electrodos_E22_mediopuño2';
eventos = [];
picos = [60];

% Baseline fija

filename = [filename '.csv'];
opts = detectImportOptions(filename);
data = readtable(filename, opts);

fs = 250; % Frecuencia de muestreo en Hz

% ---------------------------------------------------------
% Parámetros que se pueden modificar libremente
% ---------------------------------------------------------

% Tiempo para alcular promedio inicial
t_calibracion = 0.5; % Segundos de calibración inicial para quitar el offset, por defecto: 3

% Filtros
% mu es alfa, va de 7 a 14, artículos muestran 8 a 13 con mejors resultados
% beta va de 14 a 30, elegimos 13 a 30 para no chocar con mu
% emg va de 5 a 600, con fs=250 solo podemos analizar hasta 125
% mi observación empírica con emg es que es mejor de 20 a 40
% orden mayor mejor filtrado pero mayor resonancia (picos abruptos se propagan)
% orden menor peor filtrado pero menor resonancia
orden_mu = 2; % orden del filtro mu, por defecto: 2
fc1_mu = 8; fc2_mu = 13; % frecuencias de corte para mu, por defecto 8 a 13

orden_beta = 2; % orden del filtro beta, por defecto: 2
fc1_beta = 13; fc2_beta = 30; % frecuencias de corte para beta, por defecto 13 a 30

orden_emg = 2; % orden del filtro emg, por defecto: 2
fc1_emg = 20; fc2_emg = 40; % frecuencias de corte para emg, por defecto 20 a 40

orden_60 = 2; % orden del filtro notch para ruido 60Hz, por defecto 2
bw_60 = 4; % ancho de banda filtro notch, por defecto 4

eliminar_ruido_mov_resp = false; % true si se quiere eliminar ruido respiración y movimiento
orden_hp = 2; % orden del filtro pasa altas, por defecto 2
fc_hp = 1.5; % ancho de banda filtro notch, por defecto 1

% Un buen valor inicial para BCI a 250Hz es entre 0.05 y 0.15
alpha = 0.05;

% Coeficientes para la función filter: H(z) = alpha / (1 - (1-alpha)z^-1)
num_alpha = alpha;
den_alpha = [1, -(1 - alpha)];

% Envolventes de potencia
% más datos: más suave pero más retraso
% menos datos: más ruidosa pero menos retraso
n_prom_mov = round(fs/1); % 250 muestras es 1 segundo, por defecto round(fs/2) -> 0.5 segundos
%n_prom_mov = 250;
t_promedio = round(n_prom_mov/fs);
n_prom_mov_emg = round(fs/2); % por defecto: round(fs/2)

t_duracion_picos = 0.1; % Tiempo que duran los picos de interferencia
num = (1/n_prom_mov)*ones(1,n_prom_mov); % Promedia n_prom_mov datos para suavizar

% Respuesta transitoria
% En el monitor en tiempo real, debería ser casi 0 pues la señal se filtra
% desde un inicio, y se empiezan a guardar datos ya filtrados
% Aquí solo es para simular en matlab
t_resp_tr = 0.5; % tiempo estimado respuesta transitoria, por defecto: 0.5

% Umbrales de disparo
% pueden variar mucho de paciente a paciente
% los estudios son valores promedio
k_erd_mu = 0.2; % por defecto caida del 20% -> 0.2
k_erd_beta = 0.15; % por defecto caida del 15% -> 0.15
k_ers_beta = 0.2; % por defecto rebote del 20% -> 0.2

k_emg = 0.20; % por defecto 20% más energía movimiento voluntario -> 0.2
k_emg_on  = 5; % Enciende si supera 5 desviaciones estándar
k_emg_off = 4; % Apaga si baja de 2 desviaciones estándar

% Para la detección de artefactos
% mientras más se acerca a 1, más parecidas son las señales
umbral_correlacion = 0.99; % Similitud con el EMG = Artefacto, por defecto 0.6
Umbral_Ruido_Extremo = log(1 + 2.0); % picos que superan el 200% 

% Ventana móvil
% Quizás lo más importante
% ventana más grande, mejor promedio si no hay artefactos presentes
% ventana más pequeña, peor promedio pero menos posibilidad de artefactos
% pasos step_n más pequeños = detecciones más rápidas pero con más ruido
% pasos step_n más grandes = detecciones más retrasadas pero con menos ruido
t_base   = 3; % ventana base (s), por defecto 7 segundos
t_gap    = 0.25; % separación base ↔ actual (s), por defecto 1 segundo
t_actual = 3; % ventana actual (s), por defecto 2 segundos
t_step = 0.25; % Desplazamiento de ventana actual por defecto: fs * 0.25 seg (63 muestras)

n_base   = fs * t_base; % 250 muestras por cada segundo
n_gap    = round(fs * t_gap); % 250 muestras por cada segundo
n_actual = round(fs * t_actual); % 250 muestras por cada segundo
step_n = round(fs * t_step); % 250 muestras por cada segundo

tiempo_inicio = t_calibracion + t_promedio + t_base + t_gap + t_actual;

K_persistencia = 2; % Número de veces seguidas ERD ocurre, elimina falsos positivos

% ---------------------------------------------------------
% 2. Extraer los canales de interés y recrear el vector de tiempo
% ---------------------------------------------------------
fs = 250; % Frecuencia de muestreo en Hz
VRef        = 4.5; % se cambió a 5.0V
Gain        = 24.0;
Resolution  = 2^23 - 1;
Gv = (VRef / (Gain * Resolution)) * 1e6;
C3np_raw = data.CH1_C3 * Gv; % hemisferio izquierdo, músculo derecho, C3
C3p_raw = data.CH2_C3p * Gv; % hemisferio izquierdo, músculo derecho, C3'
C3pp_raw = data.CH7_C3pp * Gv; % hemisferio izquierdo, músculo derecho, C3''
C4np_raw = data.CH3_C4 * Gv; % hemisferio derecho, músculo izquierdo, C4
C4p_raw = data.CH4_C4p * Gv; % hemisferio derecho, músculo izquierdo, C4'
C4pp_raw = data.CH8_C4pp * Gv; % hemisferio derecho, músculo izquierdo, C4''
M1_raw = data.CH5_M1 * Gv; % Músculo 1
M2_raw = data.CH6_M2 * Gv; % Músculo 2

M_dif = M1_raw - M2_raw; % EMG diferencial raw
C3_raw = C3np_raw - (C3p_raw + C3pp_raw) / 2; % EEG diferencial C3 raw
C4_raw = C4np_raw - (C4p_raw + C4pp_raw) / 2; % EEG diferencial C4 raw
C43_raw = C4_raw - C3_raw; % EEG laterización

N_total = length(C4_raw); % Número de muestras totales
t_completo = (0:N_total-1) / fs; % Vector de tiempo completo

% ---------------------------------------------------------
% 3. Configuración de la ventana de calibración
% ---------------------------------------------------------
%t_calibracion = 3; % Segundos de calibración
t_calibracion = t_calibracion; % Segundos de calibración
n_calibracion = t_calibracion * fs; % 250 muestras = 1 segundo exacto
idx_inicio = n_calibracion + 1; % Empezamos a filtrar en la siguiente muestra

% Vector de tiempo ajustado para las gráficas (arranca en el segundo 3)
t_filt = t_completo(idx_inicio:end);

% Calcular y restar el promedio de calibración (Offset inicial)
C3_lista = C3_raw(idx_inicio:end) - mean(C3_raw(1:n_calibracion));
C4_lista = C4_raw(idx_inicio:end) - mean(C4_raw(1:n_calibracion));
C43_lista = C43_raw(idx_inicio:end) - mean(C43_raw(1:n_calibracion));
M_lista = M_dif(idx_inicio:end) - mean(M_dif(1:n_calibracion));

% ---------------------------------------------------------
% 4. Filtrado IIR Butterworth 
% EEG mu    (Banda 8-13 Hz) 
% EEG beta  (Banda 13-30 Hz)
% EMG       (Banda 20-40 Hz)
% ---------------------------------------------------------

% Podemos eliminar picos de frecuencias de interferencia (artefactos)

if eliminar_ruido_mov_resp
    [z, p, k] = butter(orden_hp, fc_hp/(fs/2), 'high');
    [SOS_X, G_X] = zp2sos(z, p, k);
    C3_lista = sosfilt(SOS_X, C3_lista) * prod(G_X);
    C4_lista = sosfilt(SOS_X, C4_lista) * prod(G_X);
    C43_lista = sosfilt(SOS_X, C43_lista) * prod(G_X);
    M_lista = sosfilt(SOS_X, M_lista) * prod(G_X);
    
    [Xejw,f]=freqz(SOS_X,1e4,fs);
    plot(f,abs(Xejw)*G_X,'LineWidth',1.5)
    box off; grid on; grid minor
    title('Filtro Pasa Altas');
    xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')
end

for i = 1: length(picos)
    fc1=1;
    fc2=1;
    if picos(i) == 60
        fc1=bw_60/2;
        fc2=bw_60/2;
    elseif picos(i) == 1
        fc1=0.5;
        fc2=0.5;
    elseif picos(i) == 9.76
        fc1=1;
        fc2=1;
    end
    % % Notch X Hz
    [z, p, k] = butter(orden_60, [picos(i)-fc1 picos(i)+fc1]/(fs/2), 'stop');
    [SOS_X, G_X] = zp2sos(z, p, k);

    [Xejw,f]=freqz(SOS_X,1e4,fs);
    figure; plot(f,abs(Xejw)*G_X,'LineWidth',1.5)
    box off; grid on; grid minor
    title('Filtro Rechaza Banda');
    xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

    C3_lista = sosfilt(SOS_X, C3_lista) * prod(G_X);
    C4_lista = sosfilt(SOS_X, C4_lista) * prod(G_X);
    C43_lista = sosfilt(SOS_X, C43_lista) * prod(G_X);
    M_lista = sosfilt(SOS_X, M_lista) * prod(G_X);
end

% EEG mu
% orden_mu = 2;
% fc1_mu = 8; fc2_mu = 13;
orden_mu = orden_mu;
fc1_mu = fc1_mu; fc2_mu = fc2_mu;
[z, p, k] = butter(orden_mu, [fc1_mu fc2_mu]/(fs/2), 'bandpass');
[SOS_mu, G_mu] = zp2sos(z, p, k);
C3_filt_mu = sosfilt(SOS_mu, C3_lista) * prod(G_mu);
C4_filt_mu = sosfilt(SOS_mu, C4_lista) * prod(G_mu);
C43_filt_mu = sosfilt(SOS_mu, C43_lista) * prod(G_mu);

[Xejw,f]=freqz(SOS_mu,1e4,fs);
figure; plot(f,abs(Xejw)*G_mu,'LineWidth',1.5)
box off; grid on; grid minor
title('Filtro Banda Mu');
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

% EEG beta
% orden_beta = 2; % orden del filtro beta
% fc1_beta = 13; fc2_mu = 30; % frecuencias de corte para beta
orden_beta = orden_beta; % orden del filtro beta
fc1_beta = fc1_beta; fc2_beta = fc2_beta; % frecuencias de corte para beta
[z, p, k] = butter(orden_beta, [fc1_beta fc2_beta]/(fs/2), 'bandpass');
[SOS_beta, G_beta] = zp2sos(z, p, k);
C3_filt_beta = sosfilt(SOS_beta, C3_lista) * prod(G_beta);
C4_filt_beta = sosfilt(SOS_beta, C4_lista) * prod(G_beta);
C43_filt_beta = sosfilt(SOS_beta, C43_lista) * prod(G_beta);

[Xejw,f]=freqz(SOS_beta,1e4,fs);
figure; plot(f,abs(Xejw)*G_beta,'LineWidth',1.5)
box off; grid on; grid minor
title('Filtro Banda Beta');
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

% EMG
% orden_emg = 2; % orden del filtro emg
% fc1_emg = 20; fc2_emg = 40; % frecuencias de corte para emg
orden_emg = orden_emg; % orden del filtro emg
fc1_emg = fc1_emg; fc2_emg = fc2_emg; % frecuencias de corte para emg
[z, p, k] = butter(orden_emg, [fc1_emg fc2_emg]/(fs/2), 'bandpass');
[SOS_emg, G_emg] = zp2sos(z, p, k);
M_filt = sosfilt(SOS_emg, M_lista) * prod(G_emg);

[Xejw,f]=freqz(SOS_emg,1e4,fs);
figure; plot(f,abs(Xejw)*G_emg,'LineWidth',1.5)
box off; grid on; grid minor
title('Filtro EMG');
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

% ---------------------------------------------------------
% 5. Cálculo de la Potencia y Envolventes mu con promedio
% ---------------------------------------------------------
% Esto cambió un poco, filtro promediador en lugar de pasabajas (ambos son pasabajas)
% Desventaja, en tiempo real se pierde medio segundo más fs/2 muestras (0.5 segundos)
%n_prom_mov = round(fs/2); % 250 muestras es 1 segundo -> 0.5 segundos
n_prom_mov = n_prom_mov; % 250 muestras es 1 segundo -> 0.5 segundos
n_prom_mov_emg = n_prom_mov_emg; % Promedia n_prom_mov_emg datos para suavizar
t_duracion_picos = t_duracion_picos; % Tiempo que duran los picos de interferencia
n_median_mov = 2 * round(t_duracion_picos * fs); % Elimina picos en ventanas de n_median_mov datos

num = num; % Promedia n_prom_mov datos para suavizar
den = 1; % Siempre es 1, no se mueve

potencia_cruda_C3_mu = C3_filt_mu.^2;
potencia_sin_picos_C3_mu = medfilt1(potencia_cruda_C3_mu, n_median_mov);
env_suave_C3_mu = filter(num_alpha, den_alpha, potencia_sin_picos_C3_mu);
env_log_C3_mu = log(env_suave_C3_mu+1e-8);
%C3_pwr_env_mu = env_log_C3_mu;
C3_pwr_env_mu = filter(num_alpha, den_alpha, env_log_C3_mu);

potencia_cruda_C4_mu = C4_filt_mu.^2;
potencia_sin_picos_C4_mu = medfilt1(potencia_cruda_C4_mu, n_median_mov);
env_suave_C4_mu = filter(num_alpha, den_alpha, potencia_sin_picos_C4_mu);
env_log_C4_mu = log(env_suave_C4_mu+1e-8);
%C4_pwr_env_mu = env_log_C4_mu;
C4_pwr_env_mu = filter(num_alpha, den_alpha, env_log_C4_mu);

potencia_cruda_C43_mu = C43_filt_mu.^2;
potencia_sin_picos_C43_mu = medfilt1(potencia_cruda_C43_mu, n_median_mov);
env_suave_C43_mu = filter(num_alpha, den_alpha, potencia_sin_picos_C43_mu);
env_log_C43_mu = log(env_suave_C43_mu+1e-8);
%C43_pwr_env_mu = env_log_C43_mu;
C43_pwr_env_mu = filter(num_alpha, den_alpha, env_log_C43_mu);

potencia_cruda_C3_beta = C3_filt_beta.^2;
potencia_sin_picos_C3_beta = medfilt1(potencia_cruda_C3_beta, n_median_mov);
env_suave_C3_beta = filter(num_alpha, den_alpha, potencia_sin_picos_C3_beta);
env_log_C3_beta = log(env_suave_C3_beta+1e-8);
%C3_pwr_env_beta = env_log_C3_beta;
C3_pwr_env_beta = filter(num_alpha, den_alpha, env_log_C3_beta);

potencia_cruda_C4_beta = C4_filt_beta.^2;
potencia_sin_picos_C4_beta = medfilt1(potencia_cruda_C4_beta, n_median_mov);
env_suave_C4_beta = filter(num_alpha, den_alpha, potencia_sin_picos_C4_beta);
env_log_C4_beta = log(env_suave_C4_beta+1e-8);
%C4_pwr_env_beta = env_log_C4_beta;
C4_pwr_env_beta = filter(num_alpha, den_alpha, env_log_C4_beta);

potencia_cruda_C43_beta = C43_filt_beta.^2;
potencia_sin_picos_C43_beta = medfilt1(potencia_cruda_C43_beta, n_median_mov);
env_suave_C43_beta = filter(num_alpha, den_alpha, potencia_sin_picos_C43_beta);
env_log_C43_beta = log(env_suave_C43_beta+1e-8);
%C43_pwr_env_beta = env_log_C43_beta;
C43_pwr_env_beta = filter(num_alpha, den_alpha, env_log_C43_beta);

potencia_cruda_M = M_filt.^2;
potencia_sin_picos_M = medfilt1(potencia_cruda_M, n_median_mov);
env_suave_M = filter(num_alpha, den_alpha, potencia_sin_picos_M);
env_log_M = log(env_suave_M+1e-8);
%M_pwr_env = env_log_M;
M_pwr_env = filter(num_alpha, den_alpha, env_log_M);

% ---------------------------------------------------------
% 6. Simulación en tiempo real (Simplificada con Ventanas)
% ---------------------------------------------------------

% Calcular Media y Desviación Estándar durante la calibración inicial

t_resp_tr = t_resp_tr; % tiempo estimado respuesta transitoria

idx_cal_fin = round(tiempo_inicio * fs - n_actual - n_gap);
idx_cal_ini = idx_cal_fin - n_base + 1;

C3_reposo_prom_mu = mean(C3_pwr_env_mu(idx_cal_ini : idx_cal_fin));
C3_reposo_std_mu  = std(C3_pwr_env_mu(idx_cal_ini : idx_cal_fin));
C4_reposo_prom_mu = mean(C4_pwr_env_mu(idx_cal_ini : idx_cal_fin));
C4_reposo_std_mu  = std(C4_pwr_env_mu(idx_cal_ini : idx_cal_fin));
C43_reposo_prom_mu = mean(C43_pwr_env_mu(idx_cal_ini : idx_cal_fin));
C43_reposo_std_mu  = std(C43_pwr_env_mu(idx_cal_ini : idx_cal_fin));

C3_reposo_prom_beta = mean(C3_pwr_env_beta(idx_cal_ini : idx_cal_fin));
C3_reposo_std_beta  = std(C3_pwr_env_beta(idx_cal_ini : idx_cal_fin));
C4_reposo_prom_beta = mean(C4_pwr_env_beta(idx_cal_ini : idx_cal_fin));
C4_reposo_std_beta  = std(C4_pwr_env_beta(idx_cal_ini : idx_cal_fin));
C43_reposo_prom_beta = mean(C43_pwr_env_beta(idx_cal_ini : idx_cal_fin));
C43_reposo_std_beta  = std(C43_pwr_env_beta(idx_cal_ini : idx_cal_fin));

% ---------------------------------------------------------
% NUEVO - baseline fija
% ---------------------------------------------------------
idx_rest = (t_filt >= 1 & t_filt < 9) | (t_filt >= 21 & t_filt < 29);
C3_reposo_prom_mu   = mean(C3_pwr_env_mu(idx_rest));
C3_reposo_prom_beta = mean(C3_pwr_env_beta(idx_rest));
C4_reposo_prom_mu   = mean(C4_pwr_env_mu(idx_rest));
C4_reposo_prom_beta = mean(C4_pwr_env_beta(idx_rest));
% ---------------------------------------------------------

M_reposo_prom = mean(M_pwr_env(idx_cal_ini : idx_cal_fin)); % 2.5 para estabilizar resp transitoria
M_reposo_std  = std(M_pwr_env(idx_cal_ini : idx_cal_fin));

% Definir % ERD y % ERS
% k_erd_mu = 0.20; % caida del 20%
% k_erd_beta = 0.15; % caida del 15%
% k_ers_beta = 0.20; % rebote del 20%
% k_emg = 0.20; % 20% más energía = movimiento voluntario
k_erd_mu = k_erd_mu; % caida del 20%
k_erd_beta = k_erd_beta; % caida del 15%
k_ers_beta = k_ers_beta; % rebote del 20%
k_emg = k_emg; % 20% más energía = movimiento voluntario

ERD_C3_mu_umbral = log(1-k_erd_mu); % signo - es caída
ERD_C4_mu_umbral = log(1-k_erd_mu);
ERD_C43_mu_umbral = log(1-k_erd_mu);
ERD_C3_beta_umbral = log(1-k_erd_beta);
ERD_C4_beta_umbral = log(1-k_erd_beta);
ERD_C43_beta_umbral = log(1-k_erd_beta);
ERS_C3_beta_umbral = log(1+k_ers_beta); % sigo + es rebote
ERS_C4_beta_umbral = log(1+k_ers_beta);
ERS_C43_beta_umbral = log(1+k_ers_beta);

% umbral_correlacion = 0.6; % Arriba del 60% de similitud con el EMG = Artefacto
umbral_correlacion = umbral_correlacion; % Arriba del 60% de similitud con el EMG = Artefacto
Umbral_Ruido_Extremo = Umbral_Ruido_Extremo; % picos que superan el 200% 

M_umbral = log(1+k_emg) + M_reposo_prom; % 20 más energía
%M_umbral = (1 + k_emg) * M_reposo_prom;
M_umbral_on  = M_reposo_prom + (k_emg_on * M_reposo_std);
M_umbral_off = M_reposo_prom + (k_emg_off * M_reposo_std);

muestras_totales = length(C3_pwr_env_mu);

% Parámetros de la ventana en MUESTRAS (fs = 250)
% win_n = fs * 1;             % Ventana de 1 segundo (250 muestras)
% tiempo_inicio = t_calibracion + 0.5 + 1; % 5 calibración + 0.5 filtro promedio + 1 ventana
% step_n = round(fs * 0.25);  % Desplazamiento de 0.25 seg (63 muestras)
%win_n = win_n;             % Ventana de 1 segundo (250 muestras)
tiempo_inicio = tiempo_inicio; % 5 calibración + 0.5 filtro promedio + 1 ventana
%step_n = step_n;  % Desplazamiento de 0.25 seg (63 muestras)

% Vectores de salida (inicializados en 0)
C3_erd_t_mu = zeros(muestras_totales, 1);
C4_erd_t_mu = zeros(muestras_totales, 1);
C43_erd_t_mu = zeros(muestras_totales, 1);
C3_mu_t  = zeros(muestras_totales, 1);
C4_mu_t  = zeros(muestras_totales, 1);
C43_mu_t  = zeros(muestras_totales, 1);
C3_frec_t_mu = zeros(muestras_totales, 1);
C4_frec_t_mu = zeros(muestras_totales, 1);
C43_frec_t_mu = zeros(muestras_totales, 1);

C3_erd_t_beta = zeros(muestras_totales, 1);
C4_erd_t_beta = zeros(muestras_totales, 1);
C43_erd_t_beta = zeros(muestras_totales, 1);
C3_beta_t  = zeros(muestras_totales, 1);
C4_beta_t  = zeros(muestras_totales, 1);
C43_beta_t  = zeros(muestras_totales, 1);
C3_frec_t_beta = zeros(muestras_totales, 1);
C4_frec_t_beta = zeros(muestras_totales, 1);
C43_frec_t_beta = zeros(muestras_totales, 1);

EMG_pulso_t = zeros(muestras_totales, 1);

C3_motor_t = zeros(muestras_totales, 1);
C4_motor_t = zeros(muestras_totales, 1);
C43_motor_t = zeros(muestras_totales, 1);

% Arreglos para guardar la historia de la línea base y umbrales en el tiempo
Base_C3_mu_t = zeros(muestras_totales, 1);
Base_C4_mu_t = zeros(muestras_totales, 1);
Base_C43_mu_t = zeros(muestras_totales, 1);

Base_C3_beta_t = zeros(muestras_totales, 1);
Base_C4_beta_t = zeros(muestras_totales, 1);
Base_C43_beta_t = zeros(muestras_totales, 1);

% Nuevos arreglos para guardar el estado de los artefactos
C3_artefacto_t = zeros(muestras_totales, 1);
C4_artefacto_t = zeros(muestras_totales, 1);
C43_artefacto_t = zeros(muestras_totales, 1);
artefacto_c3 = 0; artefacto_c4 = 0; artefacto_c43 = 0;

% Variables retenedoras (Zero-order hold)
erd_c3_actual_mu = 0; erd_c4_actual_mu = 0; erd_c43_actual_mu = 0;
est_c3_actual_mu = 0; est_c4_actual_mu = 0; est_c43_actual_mu = 0;
frec_mov_c3_mu = 0; frec_mov_c4_mu = 0; frec_mov_c43_mu = 0;

erd_c3_actual_beta = 0; erd_c4_actual_beta = 0; erd_c43_actual_beta = 0;
est_c3_actual_beta = 0; est_c4_actual_beta = 0; est_c43_actual_beta = 0;
frec_mov_c3_beta = 0; frec_mov_c4_beta = 0; frec_mov_c43_beta = 0;

estado_mov_emg = 0; % 0 = reposo, 1 = movimiento activo

comando_motor_c3 = 0; comando_motor_c4 = 0; comando_motor_c43 = 0;
contador_c3 = 0; contador_c4 = 0; contador_c43 = 0;

timer_c3 = 0; timer_c4 = 0; 
max_ciclos_on = round(1.0 / t_step); % 1.0 seg / 0.25 seg = 4 ciclos máximos

% Temporizadores de bloqueo para proteger la línea base
bloqueo_base_c3 = 0; bloqueo_base_c4 = 0; bloqueo_base_c43 = 0;
% Tiempo = Ventana Actual + Gap + Ventana Base
muestras_bloqueo = n_actual + n_gap + n_base;

% Arreglos para el semáforo del Frontend
C3_led_rojo_t = zeros(muestras_totales, 1);
C4_led_rojo_t = zeros(muestras_totales, 1);

% Variables de estado actual
led_rojo_c3 = 0; 
led_rojo_c4 = 0;

flanco_subida = 0;
ev=[];

% El bucle inicia a partir de que se llena la primera ventana de 1 segundo
%for i = win_n + tiempo_inicio * fs : muestras_totales
for i = round(tiempo_inicio * fs) : muestras_totales
    
    % ¿Es momento de actualizar los cálculos? (Se ejecuta solo cada ventana actual)
    if mod(i, step_n) == 0
    %if mod(i, n_actual) == 0

        rango_actual = (i - n_actual + 1) : i;
        fin_base   = i - n_actual - n_gap;
        rango_base = (fin_base - n_base + 1) : fin_base;
        
        % 1. Extraer la ventana actual de 1 segundo (Nota: sin el doble log)
        prom_mov_c3_mu   = mean(C3_pwr_env_mu(rango_actual));
        prom_mov_c4_mu   = mean(C4_pwr_env_mu(rango_actual));
        prom_mov_c43_mu  = mean(C43_pwr_env_mu(rango_actual));
        prom_mov_c3_beta = mean(C3_pwr_env_beta(rango_actual));
        prom_mov_c4_beta = mean(C4_pwr_env_beta(rango_actual));
        prom_mov_c43_beta= mean(C43_pwr_env_beta(rango_actual));
        
        % Extraer ventana del EMG para la correlación
        %win_m = M_pwr_env(rango);
        win_m = M_pwr_env(rango_actual);
        %win_m = env_log_M(rango_actual);
        prom_mov_m = mean(win_m); % Promedio de la ventana actual del EMG

        if prom_mov_m > M_umbral_on
            estado_mov_emg = 1; % Sube el pulso
            if flanco_subida == 0
                flanco_subida = 1;
                t_emg_on = i / fs;
            end
        elseif prom_mov_m < M_umbral_off
            estado_mov_emg = 0; % Baja el pulso solo si el músculo se relajó de verdad
            if flanco_subida == 1
                flanco_subida = 0;
                t_emg_off = i / fs;
                ev = [ev; [t_emg_on t_emg_off]];
            end
        end

        % ---------------------------------------------------------
        % 6.1 DETECCIÓN DE ARTEFACTOS POR CORRELACIÓN (Pearson)
        % ---------------------------------------------------------
        win_m_z = win_m - mean(win_m);
        
        % Actualizar temporizadores de bloqueo (restando el step_n en cada ciclo)
        bloqueo_base_c3 = max(0, bloqueo_base_c3 - step_n);
        bloqueo_base_c4 = max(0, bloqueo_base_c4 - step_n);
        bloqueo_base_c43 = max(0, bloqueo_base_c43 - step_n);

        % Para C3
        win_c3_z = env_log_C3_mu(rango_actual) - mean(env_log_C3_mu(rango_actual));
        corr_c3 = sum(win_c3_z .* win_m_z) / (sqrt(sum(win_c3_z.^2)) * sqrt(sum(win_m_z.^2)));
        if isnan(corr_c3), corr_c3 = 0; end
        artefacto_c3 = abs(corr_c3) > umbral_correlacion;
        
        % Si hay artefacto HOY, congelamos la base hacia el FUTURO
        if artefacto_c3 == 1
            bloqueo_base_c3 = muestras_bloqueo; 
        end

        % Para C4
        win_c4_z = env_log_C4_mu(rango_actual) - mean(env_log_C4_mu(rango_actual));
        corr_c4 = sum(win_c4_z .* win_m_z) / (sqrt(sum(win_c4_z.^2)) * sqrt(sum(win_m_z.^2)));
        if isnan(corr_c4), corr_c4 = 0; end
        artefacto_c4 = abs(corr_c4) > umbral_correlacion;
        
        if artefacto_c4 == 1
            bloqueo_base_c4 = muestras_bloqueo;
        end

        % Para C4 - C3
        win_c43_z = env_log_C43_mu(rango_actual) - mean(env_log_C43_mu(rango_actual));
        corr_c43 = sum(win_c43_z .* win_m_z) / (sqrt(sum(win_c43_z.^2)) * sqrt(sum(win_m_z.^2)));
        if isnan(corr_c43), corr_c43 = 0; end
        artefacto_c43 = abs(corr_c43) > umbral_correlacion;
        
        if artefacto_c43 == 1
            bloqueo_base_c43 = muestras_bloqueo;
        end

        % ---------------------------------------------------------
        % 6.2 CALCULAR ERD (POTENCIA LOGARÍTMICA)
        % ---------------------------------------------------------
        % Al ser logaritmos, la resta directa indica la proporción de cambio
        erd_c3_actual_mu = prom_mov_c3_mu - C3_reposo_prom_mu;
        erd_c4_actual_mu = prom_mov_c4_mu - C4_reposo_prom_mu;
        erd_c43_actual_mu = prom_mov_c43_mu - C43_reposo_prom_mu;

        erd_c3_actual_beta = prom_mov_c3_beta - C3_reposo_prom_beta;
        erd_c4_actual_beta = prom_mov_c4_beta - C4_reposo_prom_beta;
        erd_c43_actual_beta = prom_mov_c43_beta - C43_reposo_prom_beta;

        % ---------------------------------------------------------
        % 6.3 DETECCIÓN DE ARTEFACTOS POR PICOS GIGANTES
        % ---------------------------------------------------------        

        % Umbral de ruido extremo (ejemplo: si la potencia sube 200% sobre la base)
        Umbral_Ruido_Extremo = Umbral_Ruido_Extremo; 
        
        % Si AMBAS bandas experimentan una explosión de energía hacia ARRIBA al mismo tiempo:
        if (erd_c3_actual_mu > Umbral_Ruido_Extremo) && (erd_c3_actual_beta > Umbral_Ruido_Extremo)
            artefacto_c3 = 1; % ¡Es un jalón de cable o un parpadeo gigante!
            bloqueo_base_c3 = muestras_bloqueo;
        end

        % Si AMBAS bandas experimentan una explosión de energía hacia ARRIBA al mismo tiempo:
        if (erd_c4_actual_mu > Umbral_Ruido_Extremo) && (erd_c4_actual_beta > Umbral_Ruido_Extremo)
            artefacto_c4 = 1; % ¡Es un jalón de cable o un parpadeo gigante!
            bloqueo_base_c4 = muestras_bloqueo;
        end

        % Si AMBAS bandas experimentan una explosión de energía hacia ARRIBA al mismo tiempo:
        if (erd_c43_actual_mu > Umbral_Ruido_Extremo) && (erd_c43_actual_beta > Umbral_Ruido_Extremo)
            artefacto_c43 = 1; % ¡Es un jalón de cable o un parpadeo gigante!
            bloqueo_base_c43 = muestras_bloqueo;
        end

        % ---------------------------------------------------------
        % SEMÁFORO PARA EL FRONTEND (Lógica UI)
        % ---------------------------------------------------------
        % Verdadero (1) = No hacer movimiento (Artefacto o Base congelada)
        % Falso (0)   = Luz verde para mover
        
        led_rojo_c3 = (artefacto_c3 == 1) || (bloqueo_base_c3 > 0);
        led_rojo_c4 = (artefacto_c4 == 1) || (bloqueo_base_c4 > 0);

        % ---------------------------------------------------------
        % 6.4 LÍNEA BASE DINÁMICA (ÚLTIMOS X SEGUNDOS)
        % ---------------------------------------------------------
                
        % Actualizar reposo C3 SOLO si: Motor apagado Y el bloqueo expiró
        if (comando_motor_c3 == 0) && (bloqueo_base_c3 == 0)
            C3_reposo_prom_mu = mean(C3_pwr_env_mu(rango_base));
            C3_reposo_prom_beta = mean(C3_pwr_env_beta(rango_base));
        end
        
        % Actualizar reposo C4 SOLO si: Motor apagado Y el bloqueo expiró
        if (comando_motor_c4 == 0) && (bloqueo_base_c4 == 0)
            C4_reposo_prom_mu = mean(C4_pwr_env_mu(rango_base));
            C4_reposo_prom_beta = mean(C4_pwr_env_beta(rango_base));
        end

        % Actualizar reposo C4-C3 SOLO si: Motor apagado Y en ESTE instante no hay artefacto
        if (comando_motor_c43 == 0) && (bloqueo_base_c43 == 0)
            C43_reposo_prom_mu = mean(C43_pwr_env_mu(rango_base));
            C43_reposo_prom_beta = mean(C43_pwr_env_beta(rango_base));
        end
        %end

        % ---------------------------------------------------------
        % NUEVO - baseline fija
        % ---------------------------------------------------------
        idx_rest = (t_filt >= 1 & t_filt < 9) | (t_filt >= 21 & t_filt < 29);
        C3_reposo_prom_mu   = mean(C3_pwr_env_mu(idx_rest));
        C3_reposo_prom_beta = mean(C3_pwr_env_beta(idx_rest));
        C4_reposo_prom_mu   = mean(C4_pwr_env_mu(idx_rest));
        C4_reposo_prom_beta = mean(C4_pwr_env_beta(idx_rest));
        % ---------------------------------------------------------
       
        % ---------------------------------------------------------
        % 6.5 MÁQUINA DE ESTADOS COMBINADA (CON TIMEOUT Y ERS)
        % ---------------------------------------------------------
        
        % Hemisferio Izquierdo (C3) -> Controla Motor Derecho
        if comando_motor_c3 == 1
            % EL MOTOR ESTÁ ENCENDIDO: Evaluamos cuándo apagarlo
            timer_c3 = timer_c3 + 1;
            
            % Apagar si: pasó 1 segundo (timeout), hay rebote ERS, o surge un artefacto
            if (timer_c3 >= max_ciclos_on) || (erd_c3_actual_beta > ERS_C3_beta_umbral) || (artefacto_c3 == 1)
                comando_motor_c3 = 0;
                contador_c3 = 0; % Reset para exigir nueva intención
                timer_c3 = 0;
            end
            
        else
            % EL MOTOR ESTÁ APAGADO: Evaluamos si debe encenderse
            if (erd_c3_actual_mu < ERD_C3_mu_umbral) && (erd_c3_actual_beta < ERD_C3_beta_umbral) && (artefacto_c3 == 0)
                contador_c3 = contador_c3 + 1;
                
                if contador_c3 >= K_persistencia
                    comando_motor_c3 = 1; % ¡Encendido!
                    timer_c3 = 0;         % Inicia el cronómetro de 1 segundo
                end
            elseif (erd_c3_actual_beta > ERS_C3_beta_umbral) || (artefacto_c3 == 1)
                contador_c3 = 0; % Resetea la intención si hay rebote o ruido
            end
        end
        
        
        % Hemisferio Derecho (C4) -> Controla Motor Izquierdo
        if comando_motor_c4 == 1
            % EL MOTOR ESTÁ ENCENDIDO: Evaluamos cuándo apagarlo
            timer_c4 = timer_c4 + 1;
            
            if (timer_c4 >= max_ciclos_on) || (erd_c4_actual_beta > ERS_C4_beta_umbral) || (artefacto_c4 == 1)
                comando_motor_c4 = 0;
                contador_c4 = 0;
                timer_c4 = 0;
            end
            
        else
            % EL MOTOR ESTÁ APAGADO: Evaluamos si debe encenderse
            if (erd_c4_actual_mu < ERD_C4_mu_umbral) && (erd_c4_actual_beta < ERD_C4_beta_umbral) && (artefacto_c4 == 0)
                contador_c4 = contador_c4 + 1;
                
                if contador_c4 >= K_persistencia
                    comando_motor_c4 = 1;
                    timer_c4 = 0;
                end
            elseif (erd_c4_actual_beta > ERS_C4_beta_umbral) || (artefacto_c4 == 1)
                contador_c4 = 0;
            end
        end
    end
    
    % Retenedor: Asignamos el valor actual a este instante de tiempo
    C3_erd_t_mu(i) = erd_c3_actual_mu;
    C4_erd_t_mu(i) = erd_c4_actual_mu;
    C43_erd_t_mu(i) = erd_c43_actual_mu;
    C3_mu_t(i)  = est_c3_actual_mu;
    C4_mu_t(i)  = est_c4_actual_mu;
    C43_mu_t(i)  = est_c43_actual_mu;

    C3_erd_t_beta(i) = erd_c3_actual_beta;
    C4_erd_t_beta(i) = erd_c4_actual_beta;
    C43_erd_t_beta(i) = erd_c43_actual_beta;
    C3_beta_t(i)  = est_c3_actual_beta;
    C4_beta_t(i)  = est_c4_actual_beta;
    C43_beta_t(i)  = est_c43_actual_beta;

    EMG_pulso_t(i) = estado_mov_emg;

    C3_motor_t(i) = comando_motor_c3;
    C4_motor_t(i) = comando_motor_c4;

    C3_artefacto_t(i) = artefacto_c3;
    C4_artefacto_t(i) = artefacto_c4;

    C3_led_rojo_t(i) = led_rojo_c3;
    C4_led_rojo_t(i) = led_rojo_c4;

    % Retenedor Zero-Order Hold para la línea base dinámica
    Base_C3_mu_t(i) = C3_reposo_prom_mu;
    Base_C4_mu_t(i) = C4_reposo_prom_mu;
    Base_C43_mu_t(i) = C43_reposo_prom_mu;
    
    Base_C3_beta_t(i) = C3_reposo_prom_beta;
    Base_C4_beta_t(i) = C4_reposo_prom_beta;
    Base_C43_beta_t(i) = C43_reposo_prom_beta;
end

% ---------------------------------------------------------
% 7. Gráficas de Resultados
% ---------------------------------------------------------
% mu
figure;
subplot(2,1,1);
plot(t_filt, C3_erd_t_mu, 'k', 'LineWidth', 1.5); hold on;
yline(ERD_C3_mu_umbral,'r--')
box off; grid on; grid minor
title('ERD% - C3 - mu');
%legend('ERD C3 %', 'Umbral','');
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

subplot(2,1,2);
plot(t_filt, C4_erd_t_mu, 'k', 'LineWidth', 1.5); hold on;
yline(ERD_C4_mu_umbral,'r--')
box off; grid on; grid minor
title('ERD% - C4 - mu');
%legend('ERD C4 %', 'Umbral','');
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

% EMG - Mu - Motor
figure;
subplot(2,1,1);
plot(t_filt, M_pwr_env, 'k', 'LineWidth', 1.5); hold on;
plot(t_filt, EMG_pulso_t * 10, 'b-', 'LineWidth', 2);
yline(M_umbral, 'r--');
box off; grid on; grid minor
title('EMG');
%legend('EMG', 'Umbral');
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

subplot(2,1,2);
plot(t_filt, C3_led_rojo_t, 'y--', 'LineWidth', 1.5); hold on;
plot(t_filt, C4_led_rojo_t, 'c--', 'LineWidth', 1.5);
plot(t_filt, C3_motor_t, 'k', 'LineWidth', 1.5);
plot(t_filt, C4_motor_t, 'b', 'LineWidth', 1.5);
box off; grid on; grid minor
title('Comando Motor - C3 y C4');
legend('Artefacto C3','Artefacto C4','C4 motor', 'C3 motor');
ylim([-0.1 1.1]);
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

% % beta
figure;
subplot(2,1,1);
plot(t_filt, C3_erd_t_beta, 'k', 'LineWidth', 1.5); hold on;
yline(ERD_C3_beta_umbral,'r--')
yline(ERS_C3_beta_umbral,'r--')
box off; grid on; grid minor
title('ERD% - C3 - beta');
%legend('ERD C3 %', 'Umbral');
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

subplot(2,1,2);
plot(t_filt, C4_erd_t_beta, 'k', 'LineWidth', 1.5); hold on;
yline(ERD_C4_beta_umbral,'r--')
yline(ERS_C4_beta_umbral,'r--')
box off; grid on; grid minor
title('ERD% - C4 - beta');
%legend('ERD C4 %', 'Umbral');
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

% C4 - C3
figure;
subplot(2,1,1);
plot(t_filt, C43_erd_t_mu, 'k', 'LineWidth', 1.5); hold on;
yline(ERD_C43_mu_umbral,'r--')
box off; grid on; grid minor
title('ERD% - C4 - C3 - mu');
%legend('ERD C3 %', 'Umbral');
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

subplot(2,1,2);
plot(t_filt, C43_erd_t_beta, 'k', 'LineWidth', 1.5); hold on;
yline(ERD_C43_beta_umbral,'r--')
yline(ERS_C43_beta_umbral,'r--')
box off; grid on; grid minor
title('ERD% - C4 - C3 - beta');
%legend('ERD C4 %', 'Umbral');
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end
% ---------------------------------------------------------
% 7. Gráficas de Resultados (Potencia Absoluta y Umbrales Dinámicos)
% ---------------------------------------------------------
% --- GRÁFICA DE ONDA MU ---
figure;
subplot(2,1,1);
plot(t_filt, C3_pwr_env_mu, 'k', 'LineWidth', 1.5); hold on; % Señal real suavizada
plot(t_filt, Base_C3_mu_t, 'b--'); % Línea base dinámica
plot(t_filt, Base_C3_mu_t + ERD_C3_mu_umbral, 'r--'); % Umbral ERD móvil
box off; grid on; grid minor
title('Potencia Logarítmica C3 - Banda Mu');
%legend('Potencia Mu', 'Línea Base Dinámica', 'Umbral de Disparo (ERD)');
ylabel('Log(Power)');
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

subplot(2,1,2);
plot(t_filt, C4_pwr_env_mu, 'k','LineWidth', 1.5); hold on;
plot(t_filt, Base_C4_mu_t, 'b--'); 
plot(t_filt, Base_C4_mu_t + ERD_C4_mu_umbral, 'r--'); 
box off; grid on; grid minor
title('Potencia Logarítmica C4 - Banda Mu');
%legend('Potencia Mu', 'Línea Base Dinámica', 'Umbral de Disparo (ERD)');
xlabel('Tiempo (s)'); ylabel('Log(Power)');
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

% --- GRÁFICA DE ONDA BETA (Con ERD y ERS) ---
figure;
subplot(2,1,1);
plot(t_filt, C3_pwr_env_beta, 'k', 'LineWidth', 1.5); hold on;
plot(t_filt, Base_C3_beta_t, 'b--'); 
plot(t_filt, Base_C3_beta_t + ERD_C3_beta_umbral, 'r--'); % Umbral ERD (Caída)
plot(t_filt, Base_C3_beta_t + ERS_C3_beta_umbral, 'r--'); % Umbral ERS (Rebote)
box off; grid on; grid minor
title('Potencia Logarítmica C3 - Banda Beta');
%legend('Potencia Beta', 'Línea Base', 'Umbral Movimiento (ERD)', 'Umbral Frenado (ERS)');
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

subplot(2,1,2);
plot(t_filt, C4_pwr_env_beta, 'k','LineWidth', 1.5); hold on;
plot(t_filt, Base_C4_beta_t, 'b--'); 
plot(t_filt, Base_C4_beta_t + ERD_C4_beta_umbral, 'r--'); 
plot(t_filt, Base_C4_beta_t + ERS_C4_beta_umbral, 'r--'); 
box off; grid on; grid minor
title('Potencia Logarítmica C4 - Banda Beta');
%legend('Potencia Beta', 'Línea Base', 'Umbral Movimiento (ERD)', 'Umbral Frenado (ERS)');
xlabel('Tiempo (s)');
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

% Gráficas de señales raw

figure; plot(t_completo,C3_raw,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C3 sin filtrar')
xlabel('tiempo, t [s]'); ylabel('C3 [\mu V]')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

[H,f] = freqz(C3_raw,1,1e5,fs);
figure; plot(f,abs(H),'k','LineWidth',1.5)
box off; grid on; grid minor
title('Respuesta en Frecuencia C3 raw')
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

figure; plot(t_completo,C4_raw,'b','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C4 sin filtrar')
xlabel('tiempo, t [s]'); ylabel('C4 [\mu V]')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

[H,f] = freqz(C4_raw,1,1e5,fs);
figure; plot(f,abs(H),'b','LineWidth',1.5)
box off; grid on; grid minor
title('Respuesta en Frecuencia C4 raw')
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

figure; plot(t_completo,M_dif,'m','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EMG sin filtrar')
xlabel('tiempo, t [s]'); ylabel('C4 [\mu V]')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

[H,f] = freqz(M_dif,1,1e5,fs);
figure; plot(f,abs(H),'m','LineWidth',1.5)
box off; grid on; grid minor
title('Respuesta en Frecuencia EMG')
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

% Gráficas de señales sin offset

figure; plot(t_filt,C3_lista,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C3 sin offset')
xlabel('tiempo, t [s]'); ylabel('C3 [\mu V]')
legend('C3 lista')

[H,f] = freqz(C3_lista,1,1e5,fs);
figure; plot(f,abs(H),'k','LineWidth',1.5)
box off; grid on; grid minor
title('Respuesta en Frecuencia C3 lista')
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

figure; plot(t_filt,C4_lista,'b','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C4 sin offset')
xlabel('tiempo, t [s]'); ylabel('C4 [\mu V]')
legend('C4 lista')

[H,f] = freqz(C4_lista,1,1e5,fs);
figure; plot(f,abs(H),'b','LineWidth',1.5)
box off; grid on; grid minor
title('Respuesta en Frecuencia C4 lista')
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

figure; plot(t_filt,M_lista,'m','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EMG sin offset')
xlabel('tiempo, t [s]'); ylabel('EMG izq [\mu V]')
legend('EMG lista')

[H,f] = freqz(M_lista,1,1e5,fs);
figure; plot(f,abs(H),'m','LineWidth',1.5)
box off; grid on; grid minor
title('Respuesta en Frecuencia EMG lista')
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

% Gráficas de señales filtradas

figure; plot(t_filt,C3_filt_mu,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C3 filtrada mu')
xlabel('tiempo, t [s]'); ylabel('C3 [\mu V]')
%legend('C3 filt mu')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

[H,f] = freqz(C3_filt_mu,1,1e5,fs);
figure; plot(f,abs(H),'k','LineWidth',1.5)
box off; grid on; grid minor
title('Respuesta en Frecuencia C3 filt mu')
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

figure; plot(t_filt,C4_filt_mu,'b','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C4 filtrada mu')
xlabel('tiempo, t [s]'); ylabel('C4 [\mu V]')
%legend('C4 filt mu')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

[H,f] = freqz(C4_filt_mu,1,1e5,fs);
figure; plot(f,abs(H),'b','LineWidth',1.5)
box off; grid on; grid minor
title('Respuesta en Frecuencia C4 filt mu')
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

figure; plot(t_filt,C3_filt_beta,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C3 filtrada beta')
xlabel('tiempo, t [s]'); ylabel('C3 [\mu V]')
%legend('C3 filt beta')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

[H,f] = freqz(C3_filt_beta,1,1e5,fs);
figure; plot(f,abs(H),'k','LineWidth',1.5)
box off; grid on; grid minor
title('Respuesta en Frecuencia C3 filt beta')
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

figure; plot(t_filt,C4_filt_beta,'b','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C4 filtrada beta')
xlabel('tiempo, t [s]'); ylabel('C4 [\mu V]')
%legend('C4 filt beta')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

[H,f] = freqz(C4_filt_beta,1,1e5,fs);
figure; plot(f,abs(H),'b','LineWidth',1.5)
box off; grid on; grid minor
title('Respuesta en Frecuencia C4 filt beta')
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

figure; plot(t_filt,M_filt,'m','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EMG izquierdo sin offset')
xlabel('tiempo, t [s]'); ylabel('EMG izq [\mu V]')
%legend('EMG filt')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

[H,f] = freqz(M_filt,1,1e5,fs);
figure; plot(f,abs(H),'m','LineWidth',1.5)
box off; grid on; grid minor
title('Respuesta en Frecuencia EMG filt')
xlabel('frecuencia, f [Hz]'); ylabel('Magnitud')

% Gráficas de potencias

% Potencias mu
colores = lines(4);
figure; hold on;
%plot(t_filt,potencia_cruda_C3_mu,'Color',colores(1,:),'LineWidth',1.5)
plot(t_filt,potencia_sin_picos_C3_mu,'Color',colores(2,:),'LineWidth',1.5)
plot(t_filt,env_suave_C3_mu,'Color',colores(3,:),'LineWidth',1.5)
plot(t_filt,env_log_C3_mu,'Color',colores(4,:),'LineWidth',1.5)
plot(t_filt,C3_pwr_env_mu,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C3 potencia mu')
xlabel('tiempo, t [s]'); ylabel('C3 [(\mu V)^2]')
%legend('C^2','mediana','prom','log','prom log')
legend('mediana','prom','log','prom log')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

figure; hold on;
%plot(t_filt,potencia_cruda_C4_mu,'Color',colores(1,:),'LineWidth',1.5)
plot(t_filt,potencia_sin_picos_C4_mu,'Color',colores(2,:),'LineWidth',1.5)
plot(t_filt,env_suave_C4_mu,'Color',colores(3,:),'LineWidth',1.5)
plot(t_filt,env_log_C4_mu,'Color',colores(4,:),'LineWidth',1.5)
plot(t_filt,C4_pwr_env_mu,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C4 potencia mu')
xlabel('tiempo, t [s]'); ylabel('C4 [(\mu V)^2]')
%legend('C^2','mediana','prom','log','prom log')
legend('mediana','prom','log','prom log')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

figure; hold on;
%plot(t_filt,potencia_cruda_C43_mu,'Color',colores(1,:),'LineWidth',1.5)
plot(t_filt,potencia_sin_picos_C43_mu,'Color',colores(2,:),'LineWidth',1.5)
plot(t_filt,env_suave_C43_mu,'Color',colores(3,:),'LineWidth',1.5)
plot(t_filt,env_log_C43_mu,'Color',colores(4,:),'LineWidth',1.5)
plot(t_filt,C43_pwr_env_mu,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C4-C3 potencia mu')
xlabel('tiempo, t [s]'); ylabel('C4-C3 [(\mu V)^2]')
%legend('C^2','mediana','prom','log','prom log')
legend('mediana','prom','log','prom log')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

% Potencias bea
figure; hold on;
%plot(t_filt,potencia_cruda_C3_beta,'Color',colores(1,:),'LineWidth',1.5)
plot(t_filt,potencia_sin_picos_C3_beta,'Color',colores(2,:),'LineWidth',1.5)
plot(t_filt,env_suave_C3_beta,'Color',colores(3,:),'LineWidth',1.5)
plot(t_filt,env_log_C3_beta,'Color',colores(4,:),'LineWidth',1.5)
plot(t_filt,C3_pwr_env_beta,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C3 potencia beta')
xlabel('tiempo, t [s]'); ylabel('C4 [(\mu V)^2]')
%legend('C^2','mediana','prom','log','prom log')
legend('mediana','prom','log','prom log')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

figure; hold on;
%plot(t_filt,potencia_cruda_C4_beta,'Color',colores(1,:),'LineWidth',1.5)
plot(t_filt,potencia_sin_picos_C4_beta,'Color',colores(2,:),'LineWidth',1.5)
plot(t_filt,env_suave_C4_beta,'Color',colores(3,:),'LineWidth',1.5)
plot(t_filt,env_log_C4_beta,'Color',colores(4,:),'LineWidth',1.5)
plot(t_filt,C4_pwr_env_beta,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C4 potencia beta')
xlabel('tiempo, t [s]'); ylabel('C4 [(\mu V)^2]')
%legend('C^2','mediana','prom','log','prom log')
legend('mediana','prom','log','prom log')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

figure; hold on;
%plot(t_filt,potencia_cruda_C43_beta,'Color',colores(1,:),'LineWidth',1.5)
plot(t_filt,potencia_sin_picos_C43_beta,'Color',colores(2,:),'LineWidth',1.5)
plot(t_filt,env_suave_C43_beta,'Color',colores(3,:),'LineWidth',1.5)
plot(t_filt,env_log_C43_beta,'Color',colores(4,:),'LineWidth',1.5)
plot(t_filt,C43_pwr_env_beta,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C4-C3 potencia beta')
xlabel('tiempo, t [s]'); ylabel('C4-C3 [(\mu V)^2]')
%legend('C^2','mediana','prom','log','prom log')
legend('mediana','prom','log','prom log')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

% Potencias EMG
figure; hold on;
%plot(t_filt,potencia_cruda_M,'Color',colores(1,:),'LineWidth',1.5)
plot(t_filt,potencia_sin_picos_M,'Color',colores(2,:),'LineWidth',1.5)
plot(t_filt,env_suave_M,'Color',colores(3,:),'LineWidth',1.5)
plot(t_filt,M_pwr_env,'m','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EMG potencia')
xlabel('tiempo, t [s]'); ylabel('EMG [(\mu V)^2]')
legend('EMG^2','mediana','prom','log')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

% Comparación potencias
figure; hold on;
%plot(t_filt,potencia_cruda_C3_beta,'Color',colores(1,:),'LineWidth',1.5)
plot(t_filt,env_log_C3_mu,'Color',colores(4,:),'LineWidth',1.5)
plot(t_filt,env_log_C3_beta,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C3 potencias mu y beta')
xlabel('tiempo, t [s]'); ylabel('C3 [(\mu V)^2]')
%legend('C^2','mediana','prom','log','prom log')
legend('C3 mu','C3 beta')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

figure; hold on;
%plot(t_filt,potencia_cruda_C3_beta,'Color',colores(1,:),'LineWidth',1.5)
plot(t_filt,env_log_C4_mu,'Color',colores(2,:),'LineWidth',1.5)
plot(t_filt,env_log_C4_beta,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C4 potencias mu y beta')
xlabel('tiempo, t [s]'); ylabel('C4 [(\mu V)^2]')
%legend('C^2','mediana','prom','log','prom log')
legend('C4 mu','C4 beta')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

figure; hold on;
%plot(t_filt,potencia_cruda_C3_beta,'Color',colores(1,:),'LineWidth',1.5)
plot(t_filt,env_log_C43_mu,'Color',colores(2,:),'LineWidth',1.5)
plot(t_filt,env_log_C43_beta,'k','LineWidth',1.5)
box off; grid on; grid minor
title('Señal EEG C4-C3 potencias mu y beta')
xlabel('tiempo, t [s]'); ylabel('C4-C3 [(\mu V)^2]')
%legend('C^2','mediana','prom','log','prom log')
legend('C4-C3 mu','C4-C3 beta')
for i = 1:size(eventos,1)
    xline(eventos(i,1),'g')
    xline(eventos(i,2),'r')
end

freq_spec = 80;

figure;
% Espectrograma de la potencia mu para ver si los picos rítmicos
% tienen una frecuencia consistente
[S, F, T] = spectrogram(C3_pwr_env_mu - mean(C3_pwr_env_mu), ...
                        round(10*fs), round(9*fs), 1024, fs);
imagesc(T, F, 10*log10(abs(S)));
axis xy; ylim([8 13]);  % muy bajas frecuencias
xlabel('tiempo (s)'); ylabel('frecuencia de modulación (Hz)');
title('Modulación de la potencia C3 mu');
colorbar

figure;
% Espectrograma de la potencia mu para ver si los picos rítmicos
% tienen una frecuencia consistente
[S, F, T] = spectrogram(C4_pwr_env_mu - mean(C4_pwr_env_mu), ...
                        round(10*fs), round(9*fs), 1024, fs);
imagesc(T, F, 10*log10(abs(S)));
axis xy; ylim([8 13]);  % muy bajas frecuencias
xlabel('tiempo (s)'); ylabel('frecuencia de modulación (Hz)');
title('Modulación de la potencia C4 mu');
colorbar

figure;
% Espectrograma de la potencia mu para ver si los picos rítmicos
% tienen una frecuencia consistente
[S, F, T] = spectrogram(C43_pwr_env_mu - mean(C43_pwr_env_mu), ...
                        round(10*fs), round(9*fs), 1024, fs);
imagesc(T, F, 10*log10(abs(S)));
axis xy; ylim([8 13]);  % muy bajas frecuencias
xlabel('tiempo (s)'); ylabel('frecuencia de modulación (Hz)');
title('Modulación de la potencia C4-C3 mu');
colorbar

figure;
% Espectrograma de la potencia mu para ver si los picos rítmicos
% tienen una frecuencia consistente
[S, F, T] = spectrogram(C3_pwr_env_beta - mean(C3_pwr_env_beta), ...
                        round(10*fs), round(9*fs), 1024, fs);
imagesc(T, F, 10*log10(abs(S)));
axis xy; ylim([13 30]);  % muy bajas frecuencias
xlabel('tiempo (s)'); ylabel('frecuencia de modulación (Hz)');
title('Modulación de la potencia C3 beta');
colorbar

figure;
% Espectrograma de la potencia mu para ver si los picos rítmicos
% tienen una frecuencia consistente
[S, F, T] = spectrogram(C4_pwr_env_beta - mean(C4_pwr_env_beta), ...
                        round(10*fs), round(9*fs), 1024, fs);
imagesc(T, F, 10*log10(abs(S)));
axis xy; ylim([13 30]);  % muy bajas frecuencias
xlabel('tiempo (s)'); ylabel('frecuencia de modulación (Hz)');
title('Modulación de la potencia C4 beta');
colorbar

figure;
% Espectrograma de la potencia mu para ver si los picos rítmicos
% tienen una frecuencia consistente
[S, F, T] = spectrogram(C43_pwr_env_beta - mean(C43_pwr_env_beta), ...
                        round(10*fs), round(9*fs), 1024, fs);
imagesc(T, F, 10*log10(abs(S)));
axis xy; ylim([13 30]);  % muy bajas frecuencias
xlabel('tiempo (s)'); ylabel('frecuencia de modulación (Hz)');
title('Modulación de la potencia C4-C3 beta');
colorbar

figure;
% Espectrograma de la potencia mu para ver si los picos rítmicos
% tienen una frecuencia consistente
[S, F, T] = spectrogram(M_pwr_env - mean(M_pwr_env), ...
                        round(10*fs), round(9*fs), 1024, fs);
imagesc(T, F, 10*log10(abs(S)));
axis xy; ylim([20 40]);  % muy bajas frecuencias
xlabel('tiempo (s)'); ylabel('frecuencia de modulación (Hz)');
title('Modulación de la potencia EMG');
colorbar