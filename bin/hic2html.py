#!/usr/bin/env python

import sys
from pathlib import Path
import os

ASSEMBLY_MODE_TEMPLATE_DATA = """
state: "1,1,1,0,0,1,NONE",
tracks: [
    {
        name: "Scaffold boundaries",
        url: `${baseURL}/BEDPE_FILE_NAME`,
    },
    {
        name: " ",
        url: `${baseURL}/BED_FILE_NAME`,
        color: "#037ffc",
    },
],
"""

if __name__ == "__main__":
    hic_file_name = os.path.basename(sys.argv[1])
    assembly_mode = sys.argv[2] == "true"

    projectDir = "/".join(__file__.split("/")[0:-1])
    html_template_path = Path(
        f"{projectDir}/report_modules/templates/hic/hic_html_template.html"
    )

    with open(html_template_path) as f:
        html_file_lines = "".join(f.readlines())

    filled_template = html_file_lines.replace("HIC_FILE_NAME", hic_file_name)

    assembly_mode_template_data = (
        ASSEMBLY_MODE_TEMPLATE_DATA.replace(
            "BEDPE_FILE_NAME",
            f"{hic_file_name.replace('.hic', '')}.bedpe",
        ).replace(
            "BED_FILE_NAME",
            f"{hic_file_name.replace('.hic', '')}.bed",
        )
        if assembly_mode
        else ""
    )

    print(
        filled_template.replace(
            "ASSEMBLY_MODE_TEMPLATE_DATA", assembly_mode_template_data
        )
    )
