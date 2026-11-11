# Robot de 6 GDL con muñeca esférica

Brazo robótico escolar de 6 grados de libertad (UPP, Ingeniería Mecatrónica), diseñado en SolidWorks 2026 y simulado en MATLAB/Simulink (Simscape Multibody).

## Motores

| Junta | Motor | Driver |
|---|---|---|
| J1, J4, J5, J6 | GIM3505-36 | GDZ3 |
| J2 | GIM6010-48 | GDZ468H |
| J3 | GIM3510-64 | GDZ34H |

## Versiones

Las versiones 1 a 3 no existen como carpeta: esos cambios se hicieron directamente sobre la v4.

| Carpeta | Qué tiene |
|---|---|
| `Robot_6GDL_v4` | Primer diseño completo: juntas tipo horquilla con chumaceras KFL y pernos de hombro. |
| `Robot_6GDL_v5` | Cambios del compañero: baleros 696ZZ y "ejefalso" en lugar de las chumaceras. |
| `Robot_6GDL_V6` | Igual que v5, con el robot de pie. |
| `Robot_6GDL_V7` | Ensamble por eslabones (base_link … link_6) y tabla D-H validada contra el CAD. |
| `Robot_6GDL_V8` | Brazo simétrico: tornillería completa, caja de electrónica (Raspberry Pi 5), STEP por eslabón y modelo en Simscape. |
| `Robot_6GDL_V9` | Relaciones de posición con el origen de cada eslabón en su marco D-H, correcciones del profesor (placas de 1/4") y Simulink siguiendo la tabla D-H. |
| `Robot_6GDL_V10` | V9 con los orígenes de los eslabones según el dibujo de la compañera; sirvió para mostrar por qué esa tabla D-H no describe el robot. Los archivos internos conservan el nombre V9_2. |
| `Robot_6GDL_V11` | **Versión actual.** Masas de motores corregidas, análisis de par de las juntas, base nueva y generador de trayectoria quíntica en Simulink. |

## Estructura de cada versión (desde V7)

- `01_Links`: subensambles de cada eslabón.
- `02_Piezas`: piezas propias, motores, baleros, tornillería y electrónica.
- `03_STEP`: exportaciones STEP.
- `04_Documentacion`: apuntes, listas y gráficas.
- `05_Matlab`: modelos de Simulink y STEP por eslabón.
- El ensamble principal está en la raíz de cada carpeta (`Robot_6GDL_Vx.SLDASM`).

## Cómo abrirlo

1. Descarga o clona el repositorio.
2. Abre el `.SLDASM` de la versión que quieras en SolidWorks 2026 o más nuevo.
3. Los modelos de Simulink (`05_Matlab`) buscan los STEP por nombre en su misma carpeta; ábrelos desde ahí con MATLAB R2025b.
