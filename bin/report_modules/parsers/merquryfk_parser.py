import base64
import os
import re

from report_modules.parsers.parsing_commons import sort_list_of_results
from tabulate import tabulate
from pathlib import Path
import pandas as pd


def load_image_as_base64_str(file_name):
    with open(file_name, "rb") as f:
        binary_fc = f.read()

    base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
    return f"data:image/png+xml;base64,{base64_utf8_str}"


def parse_merquryfk_folder(folder_name="merquryfk_outputs"):
    dir = os.getcwdb().decode()
    merqury_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(merqury_folder_path):
        return {}

    data = {"MERQURYFK": []}

    completeness_stats_paths = [
        item for item in merqury_folder_path.glob("*.completeness.stats")
    ]

    for completeness_stats_path in completeness_stats_paths:

        data_id = os.path.basename(str(completeness_stats_path)).split(
            ".completeness.stats"
        )[0]
        individual_id = data_id.split("-reads-")[0]
        haplotypes = individual_id.split("-and-")

        completeness_stats_table = pd.read_csv(completeness_stats_path, sep="\t")
        qv_stats_table = pd.read_csv(f"{folder_name}/{data_id}.qv", sep="\t")

        data["MERQURYFK"].append(
            {
                "individual_id": individual_id,
                "completeness_stats_table": completeness_stats_table.to_dict("records"),
                "completeness_stats_table_html": tabulate(
                    completeness_stats_table,
                    headers=["Assembly", "Region", "Found", "Total", "% Covered"],
                    tablefmt="html",
                    numalign="left",
                    showindex=False,
                ),
                "qv_stats_table": qv_stats_table.to_dict("records"),
                "qv_stats_table_html": tabulate(
                    qv_stats_table,
                    headers=["Assembly", "No Support", "Total", "Error %", "QV"],
                    tablefmt="html",
                    numalign="left",
                    showindex=False,
                ),
                "plots": [
                    {
                        "hap": hap,
                        "plot": load_image_as_base64_str(
                            f"{folder_name}/{data_id}.{hap}.spectra-cn.fl.png"
                        ),
                    }
                    for hap in haplotypes
                ],
            }
        )

    return {"MERQURYFK": sort_list_of_results(data["MERQURYFK"], "individual_id")}
