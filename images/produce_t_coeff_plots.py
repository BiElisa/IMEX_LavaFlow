import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def compute_t_coeff(time, T_period, T_active, T_ramp, T_shift):
    t_rem = (time + T_shift) % T_period

    if T_ramp == 0.0:
        return np.where(t_rem <= T_active, 1.0, 0.0)
    t_coeff = np.zeros_like(time)

    for i in range(len(time)):
        if t_rem[i] <= T_ramp:
            t_coeff[i] = t_rem[i] / T_ramp
        elif t_rem[i] <= T_active - T_ramp:
            t_coeff[i] = 1.0
        elif t_rem[i] <= T_active:
            t_coeff[i] = 1.0 - (t_rem[i] - T_active + T_ramp) / T_ramp
        else:
            t_coeff[i] = 0.0
    return t_coeff

time = np.linspace(0, 150, 1000)

params_list = [
    [50.0, 50.0, 0.0, 50.0],
    [50.0, 30.0, 10.0, 0.0],
    [60.0, 40.0, 5.0, 10.0],
    [50.0, 50.0, 50.0, 0.0],
    [50.0, 50.0, 25.0, 0.0],
    [50.0, 30.0, 0.0, 10.0]
]

# --- IMPOSTAZIONI GRAFICHE GLOBALI ---
plt.rcParams.update({
    'axes.titlesize': 18,      # titolo grafico
    'axes.labelsize': 16,      # nomi assi
    'xtick.labelsize': 14,     # numeri asse x
    'ytick.labelsize': 14,     # numeri asse y
    'legend.fontsize': 13,     # legenda
})


for i, params in enumerate(params_list, 1):
    fig, ax = plt.subplots(figsize=(10, 5))
    t_coeff = compute_t_coeff(time, *params)
    ax.plot(time, t_coeff, label=f'TIME_PARAM={params}', color='orange', linewidth=2)

    # --- Titoli e labels ---
    ax.set_xlabel('Time')
    ax.set_ylabel('t_coeff')
    ax.set_title(f'Case {i}')

    # --- legenda ---
    ax.legend(loc='best')#, frameon=False)

    # --- limiti assi ---
    ax.set_xlim(-5, 155)
    ax.set_ylim(-0.05, 1.05)

    # --- griglia tratteggiata e leggera ---
    ax.grid(True, linestyle='--', linewidth=0.6, color='gray', alpha=0.4)  

    # --- Rimuove i bordi del box ---
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # --- Mantiene solo gli assi x e y visibili e neri ---
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.tight_layout()

    plt.savefig(Path.cwd() / f"t_coeff_Case_{i}.png")

    plt.pause(1.3)
    plt.close(fig)
    #plt.show()
