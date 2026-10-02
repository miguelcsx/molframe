import tempfile
import unittest
from pathlib import Path

import molframe
from molframe import motif

DATA = Path(__file__).resolve().parents[2] / "crates" / "molframe-bench" / "data"

SPECIFICATION = """
[components.site]
role = "residue"
component_ids = ["SER"]
required_atoms = ["OG"]

[[constraints]]
kind = "distance"
name = "og_cb"
first = "site.OG"
second = "site.CB"
target = 1.43
tolerance = 0.2

[profile]
id = "example-1.0"
missing = "indeterminate"

[[profile.rules]]
metric = "og_cb"
comparison = "between"
values = [1.2, 1.7]
"""


def write(text, suffix=".toml"):
    handle = tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False)
    handle.write(text)
    handle.close()
    return Path(handle.name)


class MotifTests(unittest.TestCase):
    def test_a_serine_motif_is_evaluated_at_every_mapping(self):
        structure = molframe.read(DATA / "1ubq.pdb")
        path = write(SPECIFICATION)
        self.addCleanup(path.unlink)

        report = motif.evaluate(structure, path)

        self.assertGreater(len(report), 1)
        self.assertTrue(report.mapping_ambiguous)
        first = report.evaluations[0]
        self.assertEqual(first.profile, "example-1.0")
        self.assertEqual(first.verdict, "pass")
        self.assertEqual(first.constraints[0].name, "og_cb")
        self.assertAlmostEqual(first.constraints[0].value, 1.43, delta=0.1)
        self.assertEqual([e.mapping_index for e in report.evaluations], list(range(len(report))))

    def test_a_bad_specification_or_limit_raises_value_error(self):
        structure = molframe.read(DATA / "1ubq.pdb")
        bad = write("not = [valid")
        self.addCleanup(bad.unlink)
        with self.assertRaises(ValueError):
            motif.evaluate(structure, bad)
        good = write(SPECIFICATION)
        self.addCleanup(good.unlink)
        with self.assertRaises(ValueError):
            motif.evaluate(structure, good, limits=(0, 4))
        with self.assertRaises(ValueError):
            motif.evaluate(structure, good.with_suffix(".yaml"))


if __name__ == "__main__":
    unittest.main()
