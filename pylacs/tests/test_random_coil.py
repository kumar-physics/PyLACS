import math
import pytest

from pylacs.src import random_coil as rcmod


def test_atoms_contains_expected():
    rc = rcmod.RandomCoil()
    atoms = rc.atoms()
    assert atoms == ['N', 'C', 'CA', 'CB', 'H', 'HA']


def test_get_value_single_model_wishart_HIS_C():
    rc = rcmod.RandomCoil()
    # Wishart table for H has C = 174.1 (ppm) in the uploaded table
    val = rc.get_value('HIS', 'C', 'wis')
    assert math.isclose(val, 174.1, rel_tol=0, abs_tol=1e-6)


def test_get_value_avg_two_models():
    rc = rcmod.RandomCoil()
    # Average of Wishart (174.1) and Wang (174.78) for HIS C
    expected = (174.1 + 174.78) / 2.0
    val = rc.get_value('H', 'C', ['wis', 'wan'])
    assert math.isclose(val, expected, rel_tol=0, abs_tol=1e-6)


def test_get_value_validates_residue():
    rc = rcmod.RandomCoil()
    with pytest.raises(ValueError):
        rc.get_value('Xx', 'C', 'wis')


def test_get_value_validates_atom():
    rc = rcmod.RandomCoil()
    with pytest.raises(ValueError):
        rc.get_value('H', 'QX', 'wis')


def test_get_average_structure():
    rc = rcmod.RandomCoil()
    avg = rc.get_average()  # across all tables
    # Should contain many residues, each with 6 atoms (N, C, CA, CB, H, HA)
    assert isinstance(avg, dict)
    # Pick a couple to sanity-check shape
    for res in ['H', 'A', 'G']:
        if res in avg:
            assert len(avg[res]) == 6
            assert all(isinstance(x, float) for x in avg[res])
            break
