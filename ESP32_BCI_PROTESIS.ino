/*
 * FIRMWARE DE SIMULACIÓN BCI: SOLO MOTORES (sin usar la ADS1299)
 *
 * Para pruebas de simulación: el EEG se reproduce desde un CSV en el PC (Go/Wails),
 * y el PC envía SOLO el comando 'M' para mover los servos de la prótesis. Este
 * firmware NO adquiere ni transmite EEG.
 *
 * ⚠️ SEGURIDAD DE LA ADS1299:
 * La ADS1299 puede seguir físicamente conectada. Este firmware NO la usa, pero la
 * deja en un estado SEGURO y la ignora:
 *   - CS (chip select) SIEMPRE en ALTO: con CS alto el DOUT de la ADS queda en alta
 *     impedancia, así que aunque esté convirtiendo no interfiere ni se lee nada.
 *   - Se envía UNA sola vez el comando SDATAC para detener la conversión continua.
 *   - NO se hace RESET, NI se escriben registros, NI se activa lead-off / inyección
 *     de corriente. No se ejecuta ninguna rutina que pueda dañar la tarjeta.
 *
 * Ramón Felipe Miranda Hernández
 * 2026
 *
 * Comandos PC -> ESP32 (2 bytes fijos [CMD][ARG]):
 *   'M',state -> MOTORES : bit0=C3, bit1=C4 (1=100° activo, 0=150° reposo)
 *   'I',0x00  -> IDLE     : ambos servos a reposo (150°) — por compatibilidad con Detener()
 *   (cualquier otro comando se ignora; NUNCA se arranca la ADS ni se transmite EEG)
 */

#include <SPI.h>

// --- PINOUT (idéntico al firmware normal) ---
const int PIN_CS   = 5;
const int PIN_DRDY = 4;
// SPI: 18 (CLK), 19 (MISO), 23 (MOSI)

// --- Servos MG90S (control de prótesis) ---
const int PIN_SERVO_C3 = 25;  // motor C3
const int PIN_SERVO_C4 = 26;  // motor C4
const int CH_SERVO_C3  = 4;   // canal LEDC (0-15; SPI no usa LEDC)
const int CH_SERVO_C4  = 5;
const int SERVO_FREQ   = 50;  // 50 Hz (periodo 20 ms) para MG90S
const int SERVO_RES    = 16;  // resolución PWM en bits
const int ANG_REPOSO   = 125; // posición inicial (señal = 0) 150
const int ANG_ACTIVO   = 100; // prótesis movida   (señal = 1)

// --- Comando ADS1299 que sí usamos (solo para dejarla quieta) ---
const byte SDATAC = 0x11; // detener conversión/lectura continua (estado seguro)

void setup() {
  Serial.begin(921600);

  // CS en ALTO antes que nada: la ADS queda con DOUT en alta impedancia.
  pinMode(PIN_CS, OUTPUT);
  digitalWrite(PIN_CS, HIGH);
  pinMode(PIN_DRDY, INPUT); // no se usa; solo se declara para no dejar el pin flotando como salida

  // SPI mínimo solo para poder mandar UN SDATAC y dejar la ADS sin conversión continua.
  SPI.begin(18, 19, 23, 5);
  SPI.beginTransaction(SPISettings(2000000, MSBFIRST, SPI_MODE1));
  delay(500); // esperar arranque de la ADS por si está alimentada

  stopADS_segura(); // único contacto con la ADS: SDATAC

  // Inicializar servos (PWM LEDC) y dejarlos en reposo (150°).
  ledcSetup(CH_SERVO_C3, SERVO_FREQ, SERVO_RES);
  ledcAttachPin(PIN_SERVO_C3, CH_SERVO_C3);
  ledcSetup(CH_SERVO_C4, SERVO_FREQ, SERVO_RES);
  ledcAttachPin(PIN_SERVO_C4, CH_SERVO_C4);
  setMotores(0x00); // ambos a 150°
}

void loop() {
  // Solo escuchamos comandos de 2 bytes. Nunca leemos ni transmitimos EEG.
  while (Serial.available() >= 2) {
    byte cmd = Serial.read();
    byte arg = Serial.read();

    switch (cmd) {
      case 'M': setMotores(arg);  break; // bit0=C3, bit1=C4
      case 'I': setMotores(0x00); break; // reposo (compatibilidad con Detener)
      default:  /* 'L','T','R' u otros: IGNORAR (no se toca la ADS) */ break;
    }
  }
}

// -----------------------------------------------------------------------
// ADS1299: dejarla en estado seguro (un solo SDATAC, sin más)
// -----------------------------------------------------------------------
void stopADS_segura() {
  digitalWrite(PIN_CS, LOW);
  SPI.transfer(SDATAC);   // detener conversión/lectura continua
  delayMicroseconds(5);
  digitalWrite(PIN_CS, HIGH);
  // A partir de aquí CS se queda en ALTO y no se vuelve a tocar la ADS.
}

// -----------------------------------------------------------------------
// Control de servos MG90S (PWM por LEDC) — idéntico al firmware normal
// -----------------------------------------------------------------------

// anguloAUs convierte un ángulo (0-180°) a ancho de pulso en microsegundos
// (500 us = 0°, 2500 us = 180°), rango típico de servos hobby.
int anguloAUs(int ang) {
  return 500 + (int)((long)ang * 2000 / 180);
}

// usToDuty convierte microsegundos de pulso a valor de duty para LEDC
// (periodo de 20000 us a 50 Hz, resolución SERVO_RES bits).
uint32_t usToDuty(int us) {
  return (uint32_t)((uint64_t)us * ((1UL << SERVO_RES) - 1) / 20000UL);
}

void setServo(int canal, bool activo) {
  int ang = activo ? ANG_ACTIVO : ANG_REPOSO;
  ledcWrite(canal, usToDuty(anguloAUs(ang)));
}

// setMotores aplica el estado: bit0 = C3, bit1 = C4 (1 = 100°, 0 = 150°).
void setMotores(byte state) {
  setServo(CH_SERVO_C3, state & 0x01);
  setServo(CH_SERVO_C4, state & 0x02);
}
