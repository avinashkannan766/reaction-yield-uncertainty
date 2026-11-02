from chemyield.cleaning import clean_reactions
from chemyield.profile import profile


def test_profile_reports_yield_components_and_group_feasibility(raw_synth):
    df, report = clean_reactions(raw_synth)
    text = profile(df, report, ["scaffold", "role:aryl_halide", "role:ligand", "role:nonexistent"], name="syn")
    assert "## Yield" in text and "exactly 0" in text
    assert "| role:aryl_halide |" in text
    assert "| role:ligand |" in text            # 4 ligands: a row exists whether or not a split can be built
    assert "role:nonexistent" in text           # unknown group spec must not crash the profile
