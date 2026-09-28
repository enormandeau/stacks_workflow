"""Run: python3 tests/test_self_neighbours.py [workflow_directory].

Exercises the CLI with Python 3 only. Files are confined to a temporary directory.
"""
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
script = root / "00-scripts/13_impute_missing_pairwise_distances.py"
with tempfile.TemporaryDirectory(prefix="stacks-self-neighbours-") as temporary:
    work = Path(temporary)
    vcf = work / "input.vcf"
    calls = ["1/1:20:0,20:99:100,50,0", "./.:.:.:.:.", "0/0:20:20,0:99:0,50,100"]
    vcf.write_text('\n'.join([
        "##fileformat=VCFv4.3",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc",
        "\t".join(["1", "1", "1_1", "A", "C", ".", "PASS", ".", "GT:DP:AD:GQ:PL", *calls]),
    ]) + '\n')
    checks = 0
    # a sorts before b. With a tied zero distance, dropping the first row
    # removes the valid donor a instead of the missing target b.
    for donor_distance in (0, .1):
        for include_self in (False, True):
            rows = [(a, b, 0 if a == b else donor_distance if {a,b} == {"a","b"} else 1)
                    for a in ("a", "b", "c") for b in ("a", "b", "c")
                    if include_self or a != b]
            for reverse in (False, True):
                checks += 1
                distances = work / f"distances_{checks}.tsv"
                ordered = list(reversed(rows)) if reverse else rows
                distances.write_text("\n".join("\t".join(map(str, row)) for row in ordered) + "\n")
                output = work / f"output_{checks}.vcf"
                result = subprocess.run([sys.executable, str(script), str(vcf), str(distances), "1", str(output)],
                                        capture_output=True, text=True)
                assert result.returncode == 0, result.stderr
                fields = output.read_text().splitlines()[-1].split("\t")[9:]
                assert fields[0] == calls[0] and fields[2] == calls[2]
                assert fields[1].split(":")[0] == "1/1", (donor_distance, include_self, reverse, fields)
                assert "(1/3)" in result.stdout
    print(f"PASS: {checks} CLI cases for self exclusion, zero-distance ties and input order")
