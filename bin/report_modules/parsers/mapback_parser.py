import base64
import logging
import os
import re
import sys
from collections import OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

from report_modules.parsers.parsing_commons import sort_list_of_results

LOG = logging.getLogger(__name__)


def parse_bed_gc(bed_file: str) -> OrderedDict[str, list[tuple[int, int, float]]]:
    """
    Parse BED file with GC content.

    Returns dict of contig -> list of (start, end, gc_content)
    """
    gc_data: OrderedDict[str, list[tuple[int, int, float]]] = OrderedDict()
    df = pd.read_csv(bed_file, sep="\t", header=0)
    for _, row in df.iterrows():
        contig: str = row.iloc[0]
        start: int = int(row.iloc[1]) + 1  # BED is 0-indexed
        end: int = int(row.iloc[2])  # BED is end-exclusive
        gc: float = float(row.iloc[4])
        gc_data.setdefault(contig, []).append((start, end, gc))
    return gc_data


def parse_wig(wig_file: str) -> OrderedDict[str, list[tuple[int, float]]]:
    """
    Parse WIG file for coverage.

    Returns dict of contig -> list of (position, coverage)
    """
    coverage: OrderedDict[str, list[tuple[int, float]]] = OrderedDict()
    with open(wig_file) as f:
        contig: str = ""
        span: int = 0
        cursor: int = 1
        for line in f:
            if line.startswith("track type="):
                continue
            if line.startswith("fixedStep"):
                contig = line.split("chrom=")[1].split()[0].strip()
                span = int(line.split("span=")[1].split()[0].strip())
                cursor = int(line.split("start=")[1].split()[0].strip())
                coverage.setdefault(contig, [])
                continue

            if line.strip():
                cov = float(line.strip())
                coverage[contig].append((cursor, cov))

                cursor += span
    return coverage


def plot_contig_profile(
    contig: str,
    coverage: list[tuple[int, float]],
    gc_content: list[tuple[int, int, float]],
    seq_num: int,
    folder_name: str,
    mapback_filter_length_bp: int,
    delete_images: bool = True,
) -> str:
    """
    Plot coverage and GC content for a contig and save as PNG.

    Returns base64 encoded image as string
    """

    fig, ax_cov = plt.subplots(figsize=(10, 4))
    divider = make_axes_locatable(ax_cov)
    ax_gc: Axes = divider.append_axes("top", size="40%", pad=0.2, sharex=ax_cov)

    # Coverage plot
    positions: list[int] = [pos for pos, _ in coverage]
    coverages: list[float] = [cov for _, cov in coverage]

    # Assuming that the coverage has been computed with a span of 1024
    cov_filter_len: int = 2 * int(mapback_filter_length_bp / 1024 / 2) + 1

    LOG.info(
        f"Apply filter of len {cov_filter_len} to coverage data from contig {contig}"
    )

    coverages_mm: list[float] = (
        pd.Series(coverages)
        .rolling(window=cov_filter_len, center=True, min_periods=1)
        .median()
        .to_list()
    )
    ax_cov.plot(positions, coverages_mm, color="blue", label="Coverage")
    ax_cov.set_ylabel("Coverage", color="blue")
    ax_cov.set_xlabel("Position (bp)")
    ax_cov.tick_params(axis="y", labelcolor="blue")
    ax_cov.spines["top"].set_visible(False)
    ax_cov.spines["right"].set_visible(False)

    # Increase number of ticks on x-axis
    ax_cov.xaxis.set_major_locator(MaxNLocator(nbins=15))

    # GC content plot above
    gc_positions: list[int] = [start for start, _, _ in gc_content]
    gc_values: list[float] = [gc * 100.0 for _, _, gc in gc_content]

    # Assuming that the GC content has been computed with a window of size 10000
    gc_filter_len = 2 * int(mapback_filter_length_bp / 10000 / 2) + 1

    LOG.info(f"Apply filter of len {gc_filter_len} to GC data from contig {contig}")

    gc_values_mm: list[float] = (
        pd.Series(gc_values)
        .rolling(window=gc_filter_len, center=True, min_periods=1)
        .median()
        .to_list()
    )

    ax_gc.plot(gc_positions, gc_values_mm, color="green", label="GC Content")
    ax_gc.axhline(y=50, color="black", linestyle=":", linewidth=1)
    ax_gc.set_ylabel("GC Content (%)", color="green")
    ax_gc.tick_params(axis="y", labelcolor="green")
    # ax_gc.set_ylim(0.0, 1.1)
    ax_gc.get_xaxis().set_visible(False)
    ax_gc.spines["top"].set_visible(False)
    ax_gc.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.text(0.99, 0.95, contig, ha="right", va="top", transform=ax_gc.transAxes)

    file_name = f"{folder_name}/{seq_num}.png"

    plt.savefig(file_name, dpi=300)
    plt.close(fig)

    with open(file_name, "rb") as f:
        binary_fc = f.read()

    base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
    if delete_images:
        os.remove(file_name)

    return f"data:image/png+xml;base64,{base64_utf8_str}"


def plot_profile(
    wig_file: str,
    bed_file: str,
    folder_name: str,
    mapback_filter_length_bp: int,
    delete_images: bool = True,
) -> tuple[
    OrderedDict[str, list[tuple[int, int, float]]],
    OrderedDict[str, list[tuple[int, float]]],
    list[str],
]:
    gc_data: OrderedDict[str, list[tuple[int, int, float]]] = parse_bed_gc(bed_file)
    coverage_data: OrderedDict[str, list[tuple[int, float]]] = parse_wig(wig_file)

    plots: list[str] = []
    for seq_num, contig in enumerate(coverage_data.keys(), 1):
        coverage: list[tuple[int, float]] = coverage_data[contig]
        gc_content: list[tuple[int, int, float]] = gc_data.get(contig, [])
        plots.append(
            plot_contig_profile(
                contig,
                coverage,
                gc_content,
                seq_num,
                folder_name,
                mapback_filter_length_bp,
                delete_images,
            )
        )

    return (gc_data, coverage_data, plots)


def parse_mapback_folder(
    mapback_filter_length_bp: int, folder_name: str = "mapback_outputs"
):
    dir = os.getcwdb().decode()
    reports_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(reports_folder_path):
        return {}

    list_of_bed_files = reports_folder_path.glob("*.bed")

    data: dict[str, list] = {"MAPBACK": []}

    for bed_path in list_of_bed_files:
        file_tag = re.findall(
            r"([\w]+).bed",
            bed_path.name,
        )[0]

        con_wig_path = f"{folder_name}/{file_tag}.cov.wig"
        gc_data, coverage_data, plots = plot_profile(
            con_wig_path, str(bed_path), folder_name, mapback_filter_length_bp
        )

        data["MAPBACK"].append(
            {
                "hap": file_tag,
                "gc_data": gc_data,
                "coverage_data": coverage_data,
                "plots": plots,
            }
        )

    return {"MAPBACK": sort_list_of_results(data["MAPBACK"], "hap")}


if __name__ == "__main__":
    cov_wig_file: str = sys.argv[1]
    bed_file: str = sys.argv[2]
    window_len: int = int(sys.argv[3])
    output_path: str = "./"
    delete_images: bool = False

    plot_profile(cov_wig_file, bed_file, output_path, window_len, delete_images)
