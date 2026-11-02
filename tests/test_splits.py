import numpy as np
import pytest

from chemyield.splits import group_keys, make_split


def test_partitions_are_disjoint_and_complete(clean_synth):
    for kind, gb in [("random", None), ("group", "role:aryl_halide"), ("group", "scaffold")]:
        s = make_split(clean_synth, kind, seed=1, test_size=0.2, group_by=gb)
        allidx = np.concatenate(list(s.values()))
        assert len(allidx) == len(set(allidx)) == len(clean_synth)


def test_group_split_never_leaks_groups(clean_synth):
    s = make_split(clean_synth, "group", seed=3, test_size=0.2, group_by="role:aryl_halide")
    keys = group_keys(clean_synth, "role:aryl_halide")
    test_groups = {keys[i] for i in s["test"]}
    other = {keys[i] for k in ("train", "val", "cal") for i in s[k]}
    assert test_groups and not (test_groups & other)


def test_split_is_seed_deterministic_and_seed_sensitive(clean_synth):
    a = make_split(clean_synth, "group", 0, group_by="role:aryl_halide")["test"]
    b = make_split(clean_synth, "group", 0, group_by="role:aryl_halide")["test"]
    c = make_split(clean_synth, "group", 1, group_by="role:aryl_halide")["test"]
    assert np.array_equal(a, b) and not np.array_equal(a, c)


def test_too_coarse_grouping_raises(clean_synth):
    # only 4 ligands -> any single ligand is 25% of the data, above a 10% test target
    with pytest.raises(ValueError):
        make_split(clean_synth, "group", 0, test_size=0.1, group_by="role:ligand")


def test_generic_scaffold_groups_are_coarser_than_scaffold(clean_synth):
    from chemyield.splits import group_keys
    fine, coarse = set(group_keys(clean_synth, "scaffold")), set(group_keys(clean_synth, "generic_scaffold"))
    assert 1 <= len(coarse) <= len(fine)


def test_cluster_groups_are_deterministic_and_usable(clean_synth):
    from chemyield.splits import group_keys, make_split
    a, b = group_keys(clean_synth, "cluster:6"), group_keys(clean_synth, "cluster:6")
    assert a == b and 2 <= len(set(a)) <= 6
    idx = make_split(clean_synth, "group", seed=0, group_by="cluster:6", test_size=0.2)
    held = {a[i] for i in idx["test"]}
    assert held.isdisjoint({a[i] for i in idx["train"]})       # whole clusters held out
