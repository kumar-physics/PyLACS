import numpy as np
import pytest

from pylacs.src import lacs as lacsmod


def test_mad_basic_and_zero_guard():
    x = np.array([1.0, 2.0, 100.0])
    m = lacsmod.mad(x)
    assert m > 0

    # Flat array -> MAD 0, function returns 1.0 to avoid divide-by-zero
    z = np.array([5.0, 5.0, 5.0])
    assert lacsmod.mad(z) == 1.0


def test_logistic_prob_monotone_and_center():
    z = np.array([0.0, 1.0, 2.0])
    p = lacsmod.logistic_prob(z, slope=6.0, hinge=1.0)
    assert p[0] < p[1] < p[2]
    assert abs(p[1] - 0.5) < 1e-6  # at hinge=1 -> ~0.5


def test_outlier_stats_simple_threshold():
    # Big residual should be flagged; small one not
    res = np.array([0.0, 10.0])
    flags, probs = lacsmod.outlier_stats(res, scale=1.0, cutoff_k=2.5)
    assert flags == [0, 1]
    assert probs[0] < probs[1]


def test_tag_to_label():
    tag = ('1', '1', '56', 'HIS')  # (Entity_ID, Assembly_ID, Comp_index_ID, Comp_ID)
    lbl = lacsmod.tag_to_label(tag)
    assert lbl == "56H"


def test_collect_and_report_bayes_from_fake_draws():
    # Build a minimal FitResult and matching alpha samples for one atom
    fr = lacsmod.FitResult(
        slope_pos=0.0, intercept_pos=1.0, fitted_pos=[], resid_pos=[], x_pos=[], y_pos=[], tags_pos=[],
        slope_neg=0.0, intercept_neg=3.0, fitted_neg=[], resid_neg=[], x_neg=[], y_neg=[], tags_neg=[]
    )
    fits = {'c': fr}

    # Fake posterior draws centered at the two intercepts
    alpha_samples = {
        'c': {
            'pos': np.array([0.9, 1.0, 1.1], dtype=float),
            'neg': np.array([2.9, 3.0, 3.1], dtype=float),
        }
    }

    rep = lacsmod.collect_and_report_bayes(fits, alpha_samples, cutoff_k=2.5)
    assert 'offsets' in rep and 'outliers' in rep and 'offsets_bayes' in rep
    off = rep['offsets']['c']
    off_b = rep['offsets_bayes']['c']
    # Mean offset should be around 2.0 (avg of 1 and 3)
    assert abs(off - 2.0) < 1e-6
    assert pytest.approx(2.0, abs=0.2) == off_b['mean']
    assert isinstance(off_b['ci95'], list) and len(off_b['ci95']) == 2
    assert off_b['sd'] is not None


def test_run_lacs_non_bayes_with_mocks(monkeypatch, tmp_path):
    # Mock read_star to avoid pynmrstar dependency
    fake_cs = {
        '1': {
            ('1', '1', '10', 'HIS'): {'CA': 56.0, 'CB': 30.0, 'C': 175.0, 'N': 118.0, 'H': 8.3},
            ('1', '1', '11', 'ALA'): {'CA': 52.0, 'CB': 19.0, 'C': 178.0, 'N': 124.0, 'H': 8.1},
        }
    }
    monkeypatch.setattr(lacsmod, 'read_star', lambda _: fake_cs)

    # Avoid plotting
    monkeypatch.setattr(lacsmod, 'maybe_plot_all', lambda *a, **k: None)

    # Bypass heavy libs by faking the fitter
    def fake_fit(method, xvals, yvals, tags):
        return lacsmod.FitResult(
            slope_pos=0.0, intercept_pos=1.0, fitted_pos=[], resid_pos=[],
            x_pos=[], y_pos=[], tags_pos=[],
            slope_neg=0.0, intercept_neg=3.0, fitted_neg=[], resid_neg=[],
            x_neg=[], y_neg=[], tags_neg=[]
        )
    monkeypatch.setattr(lacsmod, '_fit_atom_by_method', fake_fit)

    out = lacsmod.run_lacs("dummy.str", method='tukey', data_id='T', outdir=tmp_path, plots=False)
    assert '1' in out
    assert 'offsets' in out['1']
    # Any known atom should be present (depending on available deltas)
    assert isinstance(out['1']['offsets'], dict)


def test_run_lacs_bayes_with_mocks(monkeypatch, tmp_path):
    # Mock read_star as above
    fake_cs = {
        '1': {
            ('1', '1', '10', 'HIS'): {'CA': 56.0, 'CB': 30.0, 'C': 175.0, 'N': 118.0, 'H': 8.3},
            ('1', '1', '11', 'ALA'): {'CA': 52.0, 'CB': 19.0, 'C': 178.0, 'N': 124.0, 'H': 8.1},
        }
    }
    monkeypatch.setattr(lacsmod, 'read_star', lambda _: fake_cs)
    monkeypatch.setattr(lacsmod, 'maybe_plot_all', lambda *a, **k: None)

    # Fake bayes fitter to avoid PyMC, returning draws
    def fake_bayes(xvals, yvals, tags):
        fr = lacsmod.FitResult(
            slope_pos=0.0, intercept_pos=1.0, fitted_pos=[], resid_pos=[],
            x_pos=[], y_pos=[], tags_pos=[],
            slope_neg=0.0, intercept_neg=3.0, fitted_neg=[], resid_neg=[],
            x_neg=[], y_neg=[], tags_neg=[]
        )
        draws = {'pos': np.array([0.9, 1.1]), 'neg': np.array([2.9, 3.1])}
        return fr, draws

    monkeypatch.setattr(lacsmod, '_fit_atom_bayes', fake_bayes)

    out = lacsmod.run_lacs("dummy.str", method='bayes', data_id='B', outdir=tmp_path, plots=False)
    assert '1' in out
    rep = out['1']
    assert 'offsets_bayes' in rep
    # Each present atom should have CI info
    for atom, stats in rep['offsets_bayes'].items():
        assert 'ci95' in stats and isinstance(stats['ci95'], list) and len(stats['ci95']) == 2
