# MESA Quasi-Star Models

This repository contains the MESA (Modules for Experiments in Stellar Astrophysics) inlists and source code required to simulate the evolution of **Quasi-Stars**—massive, hydrostatic envelopes powered by the accretion luminosity of a central black hole seed.

These models are based on the physics described in **Hassan et al. (2025)**.

## 🔭 Physics Context
A Quasi-Star is a hypothetical object from the early Universe ($z \approx 10-20$) formed when a massive gas cloud collapses directly into a black hole seed (or captures one). This results in a supermassive black hole seed ($\sim 100\,M_\odot$) embedded within a massive ($10^4 - 10^6\,M_\odot$) radiation-pressure supported envelope.

Unlike normal stars, the luminosity is generated not by nuclear fusion, but by the accretion of the envelope material onto the central black hole.

**Key References:**
* **Hassan et al. (2025)** - *Primary reference for these models.* [https://ui.adsabs.harvard.edu/abs/2025arXiv251018301H/abstract]
* Coughlin and Begelman (2014) -  [https://ui.adsabs.harvard.edu/abs/2024ApJ...970..158C/abstract]
* Begelma, Rossi, Armitage (2008) - [https://ui.adsabs.harvard.edu/abs/2008MNRAS.387.1649B/abstract]
* Bellinger et al. (2023) - [https://ui.adsabs.harvard.edu/abs/2023ApJ...959..113B/abstract]

---

## ⚙️ Prerequisites & Installation

### MESA Version
This code requires **MESA Revision 24.08.1**.
Running with older or significantly newer versions may cause convergence issues due to changes in the equation of state (EOS) or opacity tables.

* Official MESA website: [https://mesastar.org/](https://mesastar.org/)
* Installation Guide: [MESA SDK & Installation](https://docs.mesastar.org/en/latest/installation.html)

### Setup
1.  Ensure the `MESA_DIR` environment variable is set.
2.  Clone this repository.
3.  Compile the code using the standard MESA clean/make sequence:
    ```bash
    ./clean
    ./mk
    ```

---

## 🚀 Running the Models

To evolve the Quasi-Star model, use the following command. Note that we use a specific run script (`rn_nomodfiles`) and a specific inlist (`inlist_evolve_header`).

```bash
./rn_nomodfiles inlist_evolve_header
```

---

## Plotting Some Quantities
(might need to install mesa_reader)

 
```import matplotlib.pyplot as plt
import mesa_reader as mr

# 1. Load the history data
# Ensure 'LOGS' is the correct path to your output directory
h = mr.MesaData('LOGS/history.data')

# 2. Setup the Plotting Environment
fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# --- Left Panel: HR Diagram ---
# MESA stores log_Teff and log_L directly
ax[0].plot(h.log_Teff, h.log_L, color='firebrick', lw=2)

# Formatting 
ax[0].invert_xaxis()
ax[0].set_xlabel(r'log(Teff / K)', fontsize=12)
ax[0].set_ylabel(r'log(L / Lsun)', fontsize=12)
ax[0].set_title('Quasi-Star HR Diagram', fontsize=14)
ax[0].grid(True, linestyle='--', alpha=0.6)

# --- Right Panel: BH Mass Evolution ---
# The black hole mass is tracked in the extra history columns.
# Depending on your run_star_extras.f, this is likely 'star_mass_black_hole' or similar.
# Standard MESA 'star_age' is in years.

# Note: Check your specific history_columns.list if 'black_hole_mass' has a different name
if 'black_hole_mass' in h.bulk_names:
    bh_mass = h.black_hole_mass
else:
    # Fallback/Example name if defined in extras
    bh_mass = h.star_mass

ax[1].plot(h.star_age, bh_mass, color='black', lw=2)

ax[1].set_xlabel('Age [years]', fontsize=12)
ax[1].set_ylabel(r'Black Hole Mass [Msun]', fontsize=12)
ax[1].set_title('Seed Growth', fontsize=14)
ax[1].grid(True, linestyle='--', alpha=0.6)

# 3. Save and Show
plt.tight_layout()
plt.savefig('qs_evolution.png', dpi=300)
print("Plot saved to qs_evolution.png")
plt.show()``` 
