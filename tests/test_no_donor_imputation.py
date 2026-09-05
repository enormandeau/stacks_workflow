"""Run: python3 tests/test_no_donor_imputation.py [workflow_directory].

Tests the CLI using temporary fixtures; requires Python 3 only.
"""
import gzip
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
script = root / "00-scripts/13_impute_missing_pairwise_distances.py"
missing = "./.:7:3,4:.:."
called = "1/1:20:0,20:99:100,50,0"
with tempfile.TemporaryDirectory(prefix="stacks-no-donor-") as temporary:
    work = Path(temporary)
    distances = work / "distances.tsv"
    # Unique zero self-distances isolate this test from the neighbour-tie issue.
    distances.write_text("a\ta\t0\na\tb\t1\na\tc\t2\nb\ta\t1\nb\tb\t0\nb\tc\t1\nc\ta\t2\nc\tb\t1\nc\tc\t0\n")
    cases = [("all_missing", [missing]*3, 1, [missing]*3, 0),
             ("empty_neighbours", [missing, called, called], 0, [missing, called, called], 0),
             ("one_imputed", [missing, called, called], 1,
              ["1/1:0:0,0:0:0,0,0", called, called], 1),
             ("partial_donor", [missing, "0/.:7:3,4:.:.", called], 1,
              [missing, "0/.:7:3,4:.:.", called], 0),
             ("dot_donor", [missing, ".:7:3,4:.:.", called], 1,
              [missing, ".:7:3,4:.:.", called], 0)]
    for compressed in (False, True):
        for name, calls, neighbours, expected, count in cases:
            extension = ".vcf.gz" if compressed else ".vcf"
            source = work / (name + extension)
            output = work / (name + ".output" + extension)
            text = '\n'.join([
                "##fileformat=VCFv4.3",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
                '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc",
                "\t".join(["1", "1", "1_1", "A", "C", ".", "PASS", ".", "GT:DP:AD:GQ:PL", *calls]),
            ]) + '\n'
            opener = gzip.open if compressed else open
            with opener(source, "wt") as handle:
                handle.write(text)
            result = subprocess.run([sys.executable, str(script), str(source), str(distances),
                                     str(neighbours), str(output)], capture_output=True, text=True)
            assert result.returncode == 0, result.stderr
            with opener(output, "rt") as handle:
                lines = handle.read().splitlines()
            assert lines[-1].split("\t")[9:] == expected, name
            assert f"({count}/3)" in result.stdout, result.stdout
            assert lines[:-1] == text.splitlines()[:-1]
print("PASS: unsupported imputations stay missing; supported calls and counts; plain/gzip input")
