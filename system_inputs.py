import pandas as pd
import numpy as np
import xlwings as xw

from tff_model import run_tff_model


def main():
    # connect to the Excel workbook and sheets
    wb = xw.Book.caller()
    sheet_inputs = wb.sheets["Inputs"]
    sheet_outputs = wb.sheets["Outputs"]

    # read input data from Excel into a DataFrame
    df = sheet_inputs.range("A1").options(
        pd.DataFrame,
        header=1,
        index=False,
        expand="table"
    ).value

    # prepare to collect all simulation results
    all_results = []

    # loop over each row in the input DataFrame to run simulations
    for run_id, row in df.iterrows():
        # initial conditions 
        y0 = [
            row["C_s0"],
            row["C_d0"],
            row["V_r0"],
            row["R_f0"]
        ]

        # parameters from Excel (override defaults)
        excel_inputs = {
            "y0": y0,
            "Q_f": row["Q_f"],
            "TMP": row["TMP"],
            "V_s": row["V_s"],
            "k_d": row["k_d"],
            "beta": row["beta"],
            "alpha": row["V_d"] / (row["V_s"] + row["V_d"])
        }

        # time settings
        t_start = row["t_start"]
        t_end   = row["t_end"]
        dt      = row["dt"]

        t_eval = np.arange(t_start, t_end + dt, dt)

        # run the model
        solution = run_tff_model(
            excel_inputs=excel_inputs,
            t_span=(t_start, t_end),
            t_eval=t_eval
        )

        # take results and store in DataFrame
        C_s, C_d, V_r, R_f = solution.y
        alpha = excel_inputs["alpha"]
        C_r = (1 - alpha) * C_s + alpha * C_d

        result_df = pd.DataFrame({
            "run_id": run_id + 1,
            "time": solution.t,
            "C_s": C_s,
            "C_d": C_d,
            "C_r": C_r,
            "V_r": V_r,
            "R_f": R_f
        })

        all_results.append(result_df)

    # combine all results into a single DataFrame
    final_output = pd.concat(all_results, ignore_index=True)
    final_output = final_output.reset_index(drop=True)

    # write results back to Excel
    sheet_outputs.clear()
    sheet_outputs.range("A1").options(index=False).value = final_output
