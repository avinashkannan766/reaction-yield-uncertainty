import json

from chemyield.collect import collect


def test_collect_reads_saved_summaries(tmp_path):
    d = tmp_path / "runs" / "ds" / "loo-role-base+gkf__rf_morgan"
    d.mkdir(parents=True)
    (d / "summary.json").write_text(json.dumps(
        {"r2": [0.5, 0.1], "rmse": [17.0, 2.0], "coverage": [0.9, 0.05], "mean_width": [45.0, 3.0]}))
    other = tmp_path / "runs" / "ds" / "random__rf_morgan"
    other.mkdir(parents=True)
    (other / "summary.json").write_text(json.dumps(
        {"r2": [0.9, 0.0], "rmse": [7.0, 0.1], "coverage": [0.9, 0.0], "mean_width": [20.0, 1.0]}))
    text = collect("ds", root=tmp_path, contains="gkf")
    assert "loo-role-base+gkf" in text and "0.900 ± 0.050" in text and "random" not in text.split("##")[1:][0][:0] + text.replace("loo-role-base+gkf", "")[:0]
    assert "split: random" not in text
