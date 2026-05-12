
# -*- coding: utf-8 -*-
"""
Gráfico del potencial electrostático a lo largo del eje Z
desde un archivo .dx
"""

import griddata
import numpy as np
import matplotlib.pyplot as plt

# ================= CONFIGURACIÓN ===================
archivo_dx = r'C:\Users\dario\OneDrive\Escritorio\3.776WTR_POTENCIAL_INVERSO.dx'
titulo_grafico = "WT Resting"
# ====================================================

KT_TO_VOLTS = 0.025692

def cargar_potencial_dx(archivo_dx):
    """Carga un archivo .dx y calcula promedio y error por plano Z"""
    dx = griddata.Grid(archivo_dx)
    data = dx.grid
    z_shape = data.shape[2]
    _, _, z_grid = np.meshgrid(*dx.edges)
    z_positions = z_grid[0, 0, :]
    dz = dx.delta[2]

    # Promedio y desviación estándar por plano Z
    valores_promedio_kT = np.zeros(z_shape)
    errores_kT = np.zeros(z_shape)
    for i in range(z_shape):
        plano = data[:, :, i]
        valores_promedio_kT[i] = plano.mean()
        errores_kT[i] = plano.std()

    # Convertir a voltios
    valores_promedio_volts = valores_promedio_kT * KT_TO_VOLTS
    errores_volts = errores_kT * KT_TO_VOLTS

    # Referencia: último plano = 0 V
    valores_promedio_volts -= valores_promedio_volts[-1]

    z_centers = z_positions[:z_shape] + dz / 2
    return z_centers, valores_promedio_volts, errores_volts


# Cargar datos
z, pot, err = cargar_potencial_dx(archivo_dx)

# ====== CONFIGURACIÓN DE FUENTE GLOBAL ======
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 24
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16

# ====== GRAFICA ======
plt.figure(figsize=(8, 12))
plt.title(titulo_grafico, pad=25, fontweight='bold', color='red')

# Potencial (desde DX)
plt.errorbar(
    pot, z, xerr=err, fmt='-o', color='blue', ecolor='lightblue',
    elinewidth=2, capsize=2, markersize=4, label='Z values 1D'
)

# Etiquetas de ejes
plt.xlabel('Potential (ΔV)', fontweight='bold')
plt.ylabel('Z axis (Å)', fontweight='bold')

# Líneas punteadas rojas en Z = -24 y 24
plt.axhline(y=-24, color='red', linestyle='--', linewidth=6)
plt.axhline(y=24, color='red', linestyle='--', linewidth=6)

# Límites del eje X
plt.xlim(-0.6, 0.6)
plt.xticks(np.arange(-0.6, 0.61, 0.2))

plt.grid(True, alpha=0.3)
plt.legend(prop={'family': 'Arial', 'weight': 'normal', 'size': 18})

plt.tight_layout()

# Guardar figura
plt.savefig(r"C:\Users\dario\OneDrive\Escritorio\GRAFICAS_4us_FINALES30062024\CAMPO_ELECTRICO_Y_POTENCIAL\WTR\Potencial_WT-R.svg", format='svg', dpi=1000)
plt.show()
