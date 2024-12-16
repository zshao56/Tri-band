import numpy as np
import matplotlib.pyplot as plt
import tidy3d as td
import tidy3d.web as web
def initialize_variables():
    global medium_WO3_Bleach, medium_WO3_Color, medium_ITO_Top, \
        medium_ITO_Bottom, medium_SPE, inf_eff, wl_start, wl_end, \
            freq_start, freq_end, freqs, freq0, freqw
    medium_WO3_Bleach = td.PoleResidue.from_file("nk/WO3-bleach-modified.json")
    medium_WO3_Color = td.PoleResidue.from_file("nk/WO3-coloration-modified.json")
    medium_ITO_Top = td.PoleResidue.from_file("nk/ITO-BaF2.json")
    medium_ITO_Bottom = td.PoleResidue.from_file("nk/ITO-glass.json")
    medium_SPE = td.PoleResidue.from_file("nk/SPE.json")
    inf_eff = 10
    wl_start = 0.4
    wl_end = 18
    freq_start = td.C_0 / wl_end
    freq_end = td.C_0 / wl_start
    freqs = np.linspace(freq_start, freq_end, 100)
    freq0 = (freq_start + freq_end) / 2
    freqw = freq_end - freq_start
initialize_variables()
# global medium_WO3_Bleach, medium_WO3_Color, medium_ITO_Top, \
#         medium_ITO_Bottom, medium_SPE, inf_eff, wl_start, wl_end, \
#             freq_start, freq_end, freqs, freq0, freqw

def make_structure(t_ITO_bottom, t_WO3, w_WO3, t_SPE, t_ITO_top, medium_WO3, folder_name):
    TRI = []
    thickness = t_SPE + t_ITO_top
    WO3_bounds = w_WO3 / 2
    ITO_bottom = td.Structure(
        geometry=td.Box.from_bounds(
            rmin=(-inf_eff, -inf_eff, -t_ITO_bottom), rmax=(inf_eff, inf_eff, 0)
        ),
        medium=medium_ITO_Bottom,
    )
    SPE = td.Structure(
        geometry=td.Box.from_bounds(
            rmin=(-inf_eff, -inf_eff, 0), rmax=(inf_eff, inf_eff, t_SPE)
        ),
        medium=medium_SPE,
    )  
    WO3 = td.Structure(
        geometry=td.Box.from_bounds(
            rmin=(-WO3_bounds, -WO3_bounds, t_SPE - t_WO3), rmax=(WO3_bounds, WO3_bounds, t_SPE)
        ),
        medium=medium_WO3,
    )
    ITO_top = td.Structure(
        geometry=td.Box.from_bounds(
            rmin=(-inf_eff, -inf_eff, t_SPE), rmax=(inf_eff, inf_eff, thickness)
        ),
        medium=medium_ITO_Top,
    )  
    TRI = [WO3]

    plane_wave = td.PlaneWave(
        source_time=td.GaussianPulse(freq0=freq0, fwidth=freqw),
        size=(td.inf, td.inf, 0),
        center=(0, 0, thickness + 0.1),
        direction="-",
        pol_angle=0,
    )
    R_monitor = td.FluxMonitor(
        center=(0, 0, thickness + 0.2),
        size=(td.inf, td.inf, 0),
        freqs=freqs,
        name="R",
    )
    T_monitor = td.FluxMonitor(
        center=(0, 0, -t_ITO_bottom - 0.1),
        size=(td.inf, td.inf, 0),
        freqs=freqs,
        name="T",
    )

    freq_1 = td.C_0 / 10
    monitor_field = td.FieldMonitor(
        center=[0, 0, 0], size=[td.inf, 0, td.inf], freqs=[freq_1], name="field"
    )
    Lz = thickness + t_ITO_bottom + 0.3
    center = thickness + 0.2 - Lz / 2  # simulation domain size in z direction
    run_time = 100 / freqw  # simulation run time

    sim = td.Simulation(
        size=(0.5, 0.5, Lz + 0.1),  # simulation domain sizes in x and y directions are set to 0
        center=(0, 0, center),
        grid_spec=td.GridSpec.auto(min_steps_per_wvl=20),
        structures=TRI,
        sources=[plane_wave],
        monitors=[R_monitor, T_monitor, monitor_field],
        run_time=run_time,
        boundary_spec=td.BoundarySpec(
            x=td.Boundary.periodic(), y=td.Boundary.periodic(), z=td.Boundary.pml()
        ),  # pml is applied in the z direction
        shutoff=1e-7,
    )  # early shutoff level is decreased to 1e-7 to increase the simulation accuracy
    return sim
# return sim

def sim_plot(sim, folder_name):

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    sim.plot(x=0, ax=ax[0])
    sim.plot_grid(x=0, ax=ax[0], lw=0.4, colors="r")
    ax[0].set_xlim(-0.6, 0.6)
    sim.plot(z=0, ax=ax[1])
    sim.plot_grid(z=0, ax=ax[1], lw=0.4, colors="r")
    ax[1].set_xlim(-0.6, 0.6)
    ax[1].set_ylim(-0.4, 0.4)
    plt.savefig(f'{folder_name}/grid_sim_structure_wwo3{w_WO3}_tSPE{t_SPE}5.png')



def plot_combined_results(sim_data_B, sim_data_C, w_WO3, t_SPE, folder_name):
    plt.figure(figsize=(10, 6))

    R_B = sim_data_B["R"].flux
    T_B = -sim_data_B["T"].flux
    A_B = 1 - R_B - T_B

    plt.plot(td.C_0 / freqs, R_B, 'r-', label='R (B)')
    plt.plot(td.C_0 / freqs, T_B, 'g-', label='T (B)')
    plt.plot(td.C_0 / freqs, A_B, 'b-', label='A (B)')

    R_C = sim_data_C["R"].flux
    T_C = -sim_data_C["T"].flux
    A_C = 1 - R_C - T_C

    plt.plot(td.C_0 / freqs, R_C, 'r--', label='R (C)')
    plt.plot(td.C_0 / freqs, T_C, 'g--', label='T (C)')
    plt.plot(td.C_0 / freqs, A_C, 'b--', label='A (C)')

    plt.xlabel("Wavelength (μm)")
    plt.ylabel("Optical Property")
    plt.xscale('log')
    plt.ylim(0, 1)
    plt.legend()
    plt.title(f"Comparison of Optical Properties\nwWO3 = {w_WO3}_tSPE = {t_SPE}")
    plt.grid(True)
    plt.savefig(f"{folder_name}/results_wwo3{w_WO3}_tSPE{t_SPE}5.png")
    plt.close() 




    output_file = f"{folder_name}/results_wwo3{w_WO3}_tSPE{t_SPE}5.txt"
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("Wavelength (μm), R (B), T (B), A (B), R (C), T (C), A (C)\n")
        for wl, r_L, t_L, a_L, r_H, t_H, a_H in zip(td.C_0 / freqs, R_B, T_B, A_B, R_C, T_C,
                                                    A_C):
            f.write(f"{wl:.6f}, {r_L:.6f}, {t_L:.6f}, \
                    {a_L:.6f}, {r_H:.6f}, {t_H:.6f}, {a_H:.6f}\n")

    print(f"Results saved to {output_file}")


t_ITO_bottom = 0.5
t_WO3 = 0.3
t_ITO_top= 0.01
w_WO3 = 0.5
t_SPE = 0.4
folder_name = './results/WO3_modified'



sims_B =  make_structure(t_ITO_bottom, t_WO3, w_WO3, t_SPE, \
                            t_ITO_top, medium_WO3_Bleach, folder_name) 
sim_plot(sims_B, folder_name)
sim_data_B = web.run(sims_B, task_name=f'B_wwo3{w_WO3}_tSPE{t_SPE}')
sims_C =  make_structure(t_ITO_bottom, \
                            t_WO3, w_WO3, t_SPE, t_ITO_top, medium_WO3_Color, folder_name) 
sim_data_C = web.run(sims_C, task_name=f'C_wwo3{w_WO3}_tSPE{t_SPE}')
plot_combined_results(sim_data_B, sim_data_C, w_WO3, t_SPE, folder_name)
