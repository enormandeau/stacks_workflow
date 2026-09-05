"""Run unchanged extraction/classification scripts on allele-swapped fixtures.

Usage: python3 tests/test_ref_alt_classification.py [workflow_directory]
Requires Python 3 and Rscript on PATH; no third-party Python modules.
All generated files are confined to a temporary directory.
"""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
scripts = root / "00-scripts"
with tempfile.TemporaryDirectory(prefix="stacks-allele-swap-") as temporary:
    work = Path(temporary)
    summaries = {}
    categories = {}
    for swapped in (False, True):
        label = "swapped" if swapped else "original"
        vcf = work / (label + ".vcf")
        lines = [
            "##fileformat=VCFv4.2",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4",
        ]
        for locus, alternate_depth in enumerate((4, 6, 8, 10, 12, 14, 16), 1):
            calls = [((0, 0), (20, 0)),
                     ((0, 1), (20-alternate_depth, alternate_depth)),
                     ((0, 1), (20-alternate_depth, alternate_depth)),
                     ((1, 1), (0, 20))]
            fields = []
            for gt, ad in calls:
                if swapped:
                    gt = tuple(1-a for a in gt)
                    ad = ad[::-1]
                fields.append(f"{gt[0]}/{gt[1]}:20:{ad[0]},{ad[1]}")
            alleles = ["C", "A"] if swapped else ["A", "C"]
            lines.append("\t".join(["1", str(locus), f"{locus}_1", *alleles,
                                     ".", "PASS", ".", "GT:DP:AD", *fields]))
        vcf.write_text("\n".join(lines) + "\n")
        summary = work / (label + ".tsv")
        subprocess.run([sys.executable, str(scripts / "08_extract_snp_duplication_info.py"),
                        str(vcf), str(summary)], check=True, capture_output=True, text=True)
        subprocess.run(["Rscript", str(scripts / "09_classify_snps.R"), str(summary)],
                       check=True, capture_output=True, text=True)
        with summary.open() as handle:
            summaries[label] = list(csv.DictReader(handle, delimiter="\t"))
        with Path(str(summary) + ".categorized").open() as handle:
            categories[label] = list(csv.DictReader(handle, delimiter="\t"))
    print("ID\toriginal_ratio\tswapped_ratio\tFis\toriginal_category\tswapped_category")
    changed = 0
    for a, b, ca, cb in zip(summaries["original"], summaries["swapped"],
                           categories["original"], categories["swapped"]):
        assert a["Fis"] == b["Fis"]
        changed += ca["Category"] != cb["Category"]
        print("\t".join([a["ID"], a["MedRatio"], b["MedRatio"], a["Fis"],
                         ca["Category"], cb["Category"]]))
    print(f"Changed classifications: {changed}/{len(categories['original'])}")
    assert changed == 0, "REF/ALT swapping changed the classification"
    expected = ["duplicated", "canonical", "canonical", "canonical",
                "canonical", "canonical", "duplicated"]
    assert [x["Category"] for x in categories["original"]] == expected
    # The raw allele ratios must remain available, not be folded in the output.
    assert float(summaries["original"][-1]["MedRatio"]) == .8
