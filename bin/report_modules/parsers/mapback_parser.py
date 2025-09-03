import base64
import logging
import os
import re
from collections import OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MaxNLocator

from report_modules.parsers.parsing_commons import sort_list_of_results

LOG = logging.getLogger(__name__)
HET_STATS_WINDOW_LEN: int = 10_000
COVERAGE_SPAN_LEN = 1024
GC_WINDOW_LEN: int = 10_000


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


def parse_het_stats(
    het_stats_file: str,
) -> OrderedDict[str, list[tuple[int, int, int, float]]]:
    """
    Parse het stats file (1-based, end inclusive).

    Returns dict of contig -> list of (start, end, het_count, mean_allele_balance)
    """
    het_data: OrderedDict[str, list[tuple[int, int, int, float]]] = OrderedDict()
    if not os.path.exists(het_stats_file):
        return het_data
    df = pd.read_csv(het_stats_file, sep="\t", header=0)
    for _, row in df.iterrows():
        contig: str = row["chrom"]
        start: int = int(row["start"])
        end: int = int(row["end"])
        het_count: int = int(row["het_count"])
        mean_allele_balance: float = float(row["mean_allele_balance"])
        het_data.setdefault(contig, []).append(
            (start, end, het_count, mean_allele_balance)
        )
    return het_data


def plot_contig_profile(
    contig: str,
    coverage: list[tuple[int, float]],
    gc_content: list[tuple[int, int, float]],
    seq_num: int,
    folder_name: str,
    mapback_rolling_median_bp: int,
    het_stats: list[tuple[int, int, int, float]] | None = None,
    delete_images: bool = True,
) -> str:
    """
    Plot coverage and GC content for a contig and save as PNG.

    Returns base64 encoded image as string
    """

    if het_stats:
        fig, axs = plt.subplots(
            4,
            1,
            figsize=(8, 6),
            sharex=True,
            gridspec_kw={"height_ratios": [1, 1, 1, 1]},
        )
        ax_gc, ax_cov, ax_het, ax_ab = axs
    else:
        fig, axs = plt.subplots(
            2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [1, 1]}
        )
        ax_gc, ax_cov = axs
        ax_het, ax_ab = None, None

    # GC content plot (top)
    gc_positions: list[int] = [start + GC_WINDOW_LEN // 2 for start, _, _ in gc_content]
    gc_values: list[float] = [gc * 100.0 for _, _, gc in gc_content]
    gc_filter_len = 2 * int(mapback_rolling_median_bp / GC_WINDOW_LEN / 2) + 1
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
    ax_gc.get_xaxis().set_visible(False)
    ax_gc.spines["top"].set_visible(False)
    ax_gc.spines["right"].set_visible(False)

    # Coverage plot
    positions: list[int] = [pos for pos, _ in coverage]
    coverages: list[float] = [cov for _, cov in coverage]
    cov_filter_len: int = 2 * int(mapback_rolling_median_bp / COVERAGE_SPAN_LEN / 2) + 1
    coverages_mm: list[float] = (
        pd.Series(coverages)
        .rolling(window=cov_filter_len, center=True, min_periods=1)
        .median()
        .to_list()
    )
    ax_cov.plot(positions, coverages_mm, color="blue", label="Coverage")
    ax_cov.set_ylabel("Coverage", color="blue")
    ax_cov.tick_params(axis="y", labelcolor="blue")
    ax_cov.spines["top"].set_visible(False)
    ax_cov.spines["right"].set_visible(False)
    ax_cov.xaxis.set_major_locator(MaxNLocator(nbins=15))
    ax_cov.get_xaxis().set_visible(False)

    # 0/1 GT Count plot
    if het_stats:
        het_positions: list[int] = [
            start + HET_STATS_WINDOW_LEN // 2 for start, _, _, _ in het_stats
        ]
        het_counts: list[int] = [count for _, _, count, _ in het_stats]

        het_counts_mm: list[float] = (
            pd.Series(het_counts)
            .rolling(window=gc_filter_len, center=True, min_periods=1)
            .median()
            .to_list()
        )

        ax_het.plot(het_positions, het_counts_mm, color="purple", label="0/1 GT Count")
        ax_het.set_ylabel("0/1 GT Count", color="purple")
        ax_het.tick_params(axis="y", labelcolor="purple")
        ax_het.spines["top"].set_visible(False)
        ax_het.spines["right"].set_visible(False)
        ax_het.get_xaxis().set_visible(False)

        # Mean allele balance plot
        ab_positions: list[int] = [
            start + HET_STATS_WINDOW_LEN // 2 for start, _, _, _ in het_stats
        ]
        ab_means: list[float] = [ab * 100.0 for _, _, _, ab in het_stats]

        ab_means_mm: list[float] = (
            pd.Series(ab_means)
            .rolling(window=gc_filter_len, center=True, min_periods=1)
            .median()
            .to_list()
        )

        ax_ab.plot(ab_positions, ab_means_mm, color="orange", label="Allele Balance")
        ax_ab.set_ylabel("Allele Balance", color="orange")
        ax_ab.tick_params(axis="y", labelcolor="orange")
        ax_ab.spines["top"].set_visible(False)
        ax_ab.spines["right"].set_visible(False)
        ax_ab.set_xlabel("Position (bp)")

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
    mapback_rolling_median_bp: int,
    het_stats_file: str,
    delete_images: bool = True,
) -> tuple[
    OrderedDict[str, list[tuple[int, int, float]]],
    OrderedDict[str, list[tuple[int, float]]],
    OrderedDict[str, list[tuple[int, int, int, float]]],
    list[str],
]:
    gc_data: OrderedDict[str, list[tuple[int, int, float]]] = parse_bed_gc(bed_file)
    coverage_data: OrderedDict[str, list[tuple[int, float]]] = parse_wig(wig_file)

    het_data: OrderedDict[str, list[tuple[int, int, int, float]]] = parse_het_stats(
        het_stats_file
    )

    plots: list[str] = []
    for seq_num, contig in enumerate(coverage_data.keys(), 1):
        coverage: list[tuple[int, float]] = coverage_data[contig]
        gc_content: list[tuple[int, int, float]] = gc_data.get(contig, [])
        het_stats: list[tuple[int, int, int, float]] | None = het_data.get(contig, None)
        plots.append(
            plot_contig_profile(
                contig,
                coverage,
                gc_content,
                seq_num,
                folder_name,
                mapback_rolling_median_bp,
                het_stats,
                delete_images,
            )
        )

    return (gc_data, coverage_data, het_data, plots)


def parse_mapback_folder(
    mapback_rolling_median_bp: int, folder_name: str = "mapback_outputs"
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
        het_stats_path = f"{folder_name}/{file_tag}.het.stats"
        gc_data, coverage_data, het_data, plots = plot_profile(
            con_wig_path,
            str(bed_path),
            folder_name,
            mapback_rolling_median_bp,
            het_stats_path,
        )

        data["MAPBACK"].append(
            {
                "hap": file_tag,
                "gc_data": gc_data,
                "coverage_data": coverage_data,
                "het_data": het_data,
                "plots": plots,
            }
        )

    return {"MAPBACK": sort_list_of_results(data["MAPBACK"], "hap")}
