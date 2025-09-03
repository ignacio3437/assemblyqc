#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict


def parse_bed(bed_path: str) -> list[tuple[str, int, int]]:
    intervals: list[tuple[str, int, int]] = []
    with open(bed_path) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            chrom, start, end = line.strip().split()[:3]
            intervals.append((chrom, int(start), int(end)))
    return intervals


def parse_vcf(vcf_path: str) -> list[tuple[str, int, float | None]]:
    het_sites: list[tuple[str, int, float | None]] = []
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            format_keys = fields[8].split(":")
            sample_values = fields[9].split(":")
            fmt: dict[str, str] = dict(zip(format_keys, sample_values))
            if fmt.get("GT") != "0/1":
                continue

            ad = fmt.get("AD")
            if ad:
                ad_ref, ad_alt = [int(x) for x in ad.split(",")]
                ad_ratio: float | None = (
                    ad_alt / (ad_ref + ad_alt) if (ad_ref + ad_alt) > 0 else None
                )
            else:
                ad_ratio = None
            het_sites.append((chrom, pos, ad_ratio))
    return het_sites


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Count het sites and mean allelic depth ratio per interval."
    )
    parser.add_argument("--vcf", required=True, help="Path to VCF file.")
    parser.add_argument("--bed", required=True, help="Path to BED file with intervals.")
    args = parser.parse_args()

    bed_path: str = args.bed
    vcf_path: str = args.vcf
    intervals: list[tuple[str, int, int]] = parse_bed(bed_path)
    het_sites: list[tuple[str, int, float | None]] = parse_vcf(vcf_path)

    # Index het sites by chrom for fast lookup
    sites_by_chrom: dict[str, list[tuple[int, float | None]]] = defaultdict(list)
    for chrom, pos, ad_ratio in het_sites:
        sites_by_chrom[chrom].append((pos, ad_ratio))

    print("chrom\tstart\tend\thet_count\tmean_allele_balance")
    for chrom, start, end in intervals:
        sites: list[float] = [
            ad
            for pos, ad in sites_by_chrom.get(chrom, [])
            if start <= pos < end and ad is not None
        ]
        het_count: int = len(sites)
        mean_ad: float = sum(sites) / het_count if het_count > 0 else 0.0
        print(
            f"{chrom}\t{start + 1}\t{end}\t{het_count}\t{mean_ad:.4f}"
        )  # 1-based, end inclusive


if __name__ == "__main__":
    main()
