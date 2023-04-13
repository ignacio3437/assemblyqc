import os
from pathlib import Path
import base64
import re


def parse_circos_folder(folder_name="circos_outputs"):
    dir = os.getcwdb().decode()
    circos_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(circos_folder_path):
        return {}

    list_of_plot_files = circos_folder_path.glob("*.svg")

    data = {"CIRCOS": []}

    for plot_path in list_of_plot_files:
        binary_fc = open(plot_path, "rb").read()
        base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
        ext = str(plot_path).split(".")[-1]
        plot_url = f"data:image/{ext}+xml;base64,{base64_utf8_str}"

        file_tokens = re.findall(
            r"([\w]+).on.([\w]+).svg",
            os.path.basename(str(plot_path)),
        )[0]

        data["CIRCOS"].append(
            {
                "tag.on.tag": f"{file_tokens[0]} : {file_tokens[1]}",
                "circos_plot": plot_url,
            }
        )

    if len(data["CIRCOS"]) < 1:
        return {}

    return data
