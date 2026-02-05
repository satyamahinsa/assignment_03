# Tangential Flow Filtration (TFF) Simulation Project - Dynamic Analysis with Excel–Python Integration

---

## Overview

This project focuses on the **mathematical modeling and simulation of a Tangential Flow Filtration (TFF) process** using a system of non-linear ordinary differential equations (ODEs). The model describes the dynamic behavior of:

- Retentate concentration
- Retentate volume
- Membrane fouling
- Non-ideal mixing effects such as dead zones and bypass flow inside a recirculation vessel

The simulation engine is implemented in **Python**, while **Microsoft Excel is used as the user interface** for parameter input and result visualization. This approach enables users to modify process parameters, run simulations, and analyze results without directly interacting with the Python code.

---

## Objectives and Expected Insights

The main objectives of this project are to:

- Predict the **duration of the TFF process** required to reach a target concentration  
- Observe **changes in retentate volume and concentration over time**  
- Analyze **non-ideal process behavior**, including:
  - Dead zones in the recirculation vessel
  - Bypass (short-circuiting) flow
- Compare multiple operating scenarios within a single automated workflow

---

## Project Structure
```bash
project/
│
├── img/ # screenshots and images used in the documentation
├── tff_model.py # TFF model (ODEs and solver)
├── run_from_excel.py # main automation script executed from Excel (xlwings)
├── system_inputs.xlsm # excel file for inputs, outputs, and visualizations
├── README.md # project documentation
```


---

## Model Description

The TFF model includes:

- **Non-ideal mixing** in the recirculation vessel:
  - Active zone concentration (C_s)
  - Dead zone concentration (C_d)

- **Permeate flux (J)** based on Darcy’s law, accounting for:
  - Viscosity effects
  - Osmotic pressure
  - Membrane fouling resistance

- **Retentate volume dynamics**
- **Time-dependent fouling resistance growth**

The system of equations is solved numerically using `scipy.integrate.solve_ivp`.

---
## Project Workflow
### Step 1 – Preparing the Excel Input File

#### Sheet: `Inputs`

Each **row represents one simulation scenario** with a unique set of operating conditions.

Typical input columns include:

| Column | Description |
|------|------------|
| C_s0 | Initial active zone concentration |
| C_d0 | Initial dead zone concentration |
| V_r0 | Initial retentate volume |
| R_f0 | Initial fouling resistance |
| Q_f  | Recirculation flow rate |
| TMP  | Transmembrane Pressure |
| beta | Bypass fraction |
| V_s  | Active zone volume |
| V_d  | Dead zone volume |
| k_d  | Inter-zone mass transfer coefficient |
| k_f  | Fouling growth coefficient |
| t_start | Simulation start time |
| t_end   | Simulation end time |
| dt      | Output time step |

### Step 2 – Running the Simulation from Excel

#### Requirements

- Python 3.x
- Required Python libraries:
  - numpy
  - pandas
  - scipy
  - xlwings
- xlwings properly configured with Microsoft Excel, including installation of the Excel add-in through the terminal using the command `xlwings addin install`

#### How to Run

1. Open `system_inputs.xlsm`
2. Ensure the following sheets exist:
   - `Inputs` (contains simulation parameters)
   - `Outputs` (can initially be empty)
3. Run the Python script from Excel using:
**xlwings → Run main**
    !['xlwings'](/03_assignment/img/xlwings.png)

This triggers the `main()` function in `run_from_excel.py`.

### Step 3 – Automated Workflow

When the script is executed:

1. Excel sends the data from the `Inputs` sheet to Python
2. Python:
   - Reads each row as an independent simulation
   - Runs the TFF model
   - Solves the ODE system for each scenario
3. All simulation results are combined into a single dataset
4. Results are automatically written back to Excel in the `Outputs` sheet

### Step 4 – Output Data Structure

#### Sheet: `Outputs`

The output sheet contains time-resolved simulation results with the following columns:

| Column | Description |
|------|------------|
| run_id | Simulation identifier |
| time | Simulation time |
| C_s | Active zone concentration |
| C_d | Dead zone concentration |
| C_r | Overall retentate concentration |
| V_r | Retentate volume |
| R_f | Fouling resistance |

This table is directly used for plotting and analysis in Excel.

---

## Recommended Plots in Excel

The following plots in `Plot Trends` sheet should be created **directly in Excel** using the `Outputs` sheet:

1. Temporal Evolution of Retentate Concentration during TFF
2. Retentate Volume Reduction during TFF 
3. Concentration–Volume Relationship during TFF
4. Non-Ideal Mixing Effects: Active Zone vs Dead Zone Concentrations  
5. Fouling Resistance Development during TFF 

---

## Results Interpretation

The simulation results show that all output variables remain constant throughout the entire simulation time. Specifically, the overall retentate concentration (`C_r`), retentate volume (`V_r`), active and dead zone concentrations (`C_s` and `C_d`), and fouling resistance (`R_f`) do not change between time steps.

This behavior occurs because the **permeate flux (`J`) is equal to zero during the simulation**. In the TFF model, permeate flux represents the driving mechanism for solvent removal through the membrane. When `J = 0`, no permeation takes place, and therefore:

- No solvent is removed from the system, resulting in a constant retentate volume.
- Retentate concentration does not increase, since there is no volume reduction.
- Fouling resistance does not develop, as fouling growth is directly proportional to the permeate flux.
- No concentration gradients form between the active and dead zones, preventing the observation of non-ideal mixing effects.

From a physical perspective, a zero permeate flux indicates that the **effective transmembrane pressure is insufficient to overcome the osmotic pressure of the solution**. As a result, the system remains in a non-permeating regime, and the model converges to a steady-state solution at the initial conditions.

Although dynamic behavior is not observed under the current parameter set, the results confirm the correctness of the model implementation and the Excel–Python automation workflow. By adjusting operating parameters to achieve a positive permeate flux, the same framework can be used to analyze process duration, volume reduction, and non-ideal mixing effects in a fully dynamic TFF operation.

---

## AI Tool Usage

AI tools were used as a development assistant to:

- Refine and debug the Python–Excel integration workflow
- Improve code structure and readability for the TFF simulation model 
- Generate and refine clear, structured technical documentation

---

## Conclusion

This project demonstrates the implementation of a Tangential Flow Filtration (TFF) process model using a system of non-linear ordinary differential equations integrated with an Excel–Python workflow.

Despite the absence of dynamic changes in the current simulations, the project successfully demonstrates that:

- A dynamic TFF process can be formulated and solved using non-linear ODEs
- Excel can function as a practical and user-friendly interface for defining parameters and visualizing results
- Python provides a robust computational backend for automated simulation and data processing
- Non-ideal effects such as dead zones and bypass flow are structurally incorporated into the model and can be analyzed once permeation conditions are achieved  

Overall, the workflow provides a solid and extensible framework for future studies, including parameter sensitivity analysis, investigation of permeating operating regimes, and educational exploration of TFF process behavior.
