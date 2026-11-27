"""Regression against real protease atom identities; run inside the RMSX image."""
import importlib.util
from pathlib import Path
import unittest

import MDAnalysis as mda
from rmsx.core import get_selection_string

WRAPPER = Path(__file__).resolve().parents[1] / 'tools/flipbook'
spec = importlib.util.spec_from_file_location('multichain', WRAPPER / 'rmsx_multichain.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class SharedSegmentRegression(unittest.TestCase):
    def test_ambiguous_topology_rejected_in_all_selection_modes(self):
        u = mda.Universe(str(WRAPPER / 'test-data/protease_shared_segment.pdb'))
        self.assertEqual(set(u.atoms.chainIDs), {'A', 'B'})
        self.assertEqual(set(u.atoms.segids), {'PROT'})
        for mode, selector in [('all', ''), ('selected', 'B'), ('selected', 'PROT')]:
            with self.subTest(mode=mode, selector=selector):
                with self.assertRaisesRegex(ValueError, 'Distinct chains share segment ID'):
                    module.resolve_chains(u, mode, selector, 'protein', get_selection_string)

    def test_original_protease_preserves_chain_selection(self):
        u = mda.Universe(str(WRAPPER / 'test-data/protease_backbone.pdb'))
        self.assertEqual(module.resolve_chains(u, 'all', '', 'protein', get_selection_string), ['A', 'B'])
        self.assertEqual(module.resolve_chains(u, 'selected', 'B', 'protein', get_selection_string), ['B'])

    def test_unsupported_modes_rejected(self):
        u = mda.Universe(str(WRAPPER / 'test-data/protease_backbone.pdb'))
        for mode in ['dna', 'rna', 'generic']:
            with self.subTest(mode=mode), self.assertRaisesRegex(ValueError, 'protein analysis only'):
                module.resolve_chains(u, 'all', '', mode, get_selection_string)

    def test_fixture_preserves_coordinates_and_overlapping_residues(self):
        a = mda.Universe(str(WRAPPER / 'test-data/protease_backbone.pdb'))
        b = mda.Universe(str(WRAPPER / 'test-data/protease_shared_segment.pdb'))
        self.assertTrue((a.atoms.positions == b.atoms.positions).all())
        for chain in ['A', 'B']:
            self.assertEqual(set(b.select_atoms('chainID ' + chain).resids), set(range(1, 100)))


if __name__ == '__main__':
    unittest.main()
