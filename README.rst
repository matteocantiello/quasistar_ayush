# MESA Quasi-Star Models

This repository contains the MESA (Modules for Experiments in Stellar Astrophysics) inlists and source code required to simulate the evolution of **Quasi-Stars**—massive, hydrostatic envelopes powered by the accretion luminosity of a central black hole seed.

These models are based on the physics described in **Hassan et al. (2025)**.

## 🔭 Physics Context
A Quasi-Star is a hypothetical object from the early Universe ($z ~ 10-20$) formed when a massive gas cloud collapses directly into a black hole seed (or captures one). This results in a supermassive black hole seed ($~ 100\,M_\odot$) embedded within a massive ($10^4 - 10^6\,M_\odot$) radiation-pressure supported envelope.

Unlike normal stars, the luminosity is generated not by nuclear fusion, but by the accretion of the envelope material onto the central black hole.

**Key References:**
* **Hassan et al. (2025)** - *Primary reference for these models.*
* Begelman, Volonteri, & Rees (2006) - *Foundational theory.*
* Ball, Tout, & Żytkow (2011) - *Previous MESA implementation.*

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
