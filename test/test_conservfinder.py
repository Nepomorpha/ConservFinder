import tempfile
import unittest
from pathlib import Path

from conservfinder import NtCounter, indices_to_ranges, process_alignments


class TestConservFinder(unittest.TestCase):
    def test_nt_counter(self):
        self.assertEqual(NtCounter(["AAAA", "AAAA", "TTTT"], 0.6), [0, 1, 2, 3])
        self.assertEqual(NtCounter(["----", "----", "AAAA"], 0.6), [])
        self.assertEqual(NtCounter(["NNNN", "NNNN", "AAAA"], 0.6), [])

    def test_ranges(self):
        self.assertEqual(indices_to_ranges(list(range(10))), [(0, 9)])
        self.assertEqual(indices_to_ranges(list(range(9))), [])

    def test_maf(self):
        maf = """##maf version=1

a
s sp1.chr1 10 5 + 100 --AAAAA
s sp2.chr2 20 7 + 100 TTAAAAA
s sp3.chr3 30 5 - 100 --AAAAA

a
s sp1.chr1 20 5 + 100 AAAAA
s sp2.chr2 30 5 + 100 AAAAA

a
s sp1.chr1 20 5 + 100 AAAAA
s sp1.chr4 40 5 + 100 AAAAA
s sp2.chr2 30 5 + 100 AAAAA
s sp3.chr3 40 5 - 100 AAAAA

"""
        with tempfile.TemporaryDirectory() as folder:
            input_file = Path(folder) / "test.maf"
            output_file = Path(folder) / "test.rbh2"
            input_file.write_text(maf)

            process_alignments(
                input_file,
                ["sp1.chr1", "sp2.chr2", "sp3.chr3"],
                "cactus",
                0.6,
                str(output_file),
                5
            )

            lines = output_file.read_text().splitlines()
            self.assertEqual(len(lines), 2)
            header = lines[0].split("\t")
            values = lines[1].split("\t")
            row = {}
            for i in range(len(header)):
                row[header[i]] = values[i]

            self.assertEqual(row["sp1_scaf"], "chr1")
            self.assertEqual(row["sp1_start"], "10")
            self.assertEqual(row["sp1_stop"], "15")
            self.assertEqual(row["sp2_start"], "22")
            self.assertEqual(row["sp2_stop"], "27")
            self.assertEqual(row["sp3_start"], "65")
            self.assertEqual(row["sp3_stop"], "70")
            self.assertEqual(row["sp3_strand"], "-")


if __name__ == "__main__":
    unittest.main()
