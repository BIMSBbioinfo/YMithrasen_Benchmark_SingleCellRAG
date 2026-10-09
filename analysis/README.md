# Analysis

```bash
python analysis/summarise_results.py   # results/api_summary.csv, results/coding_agents_summary.csv
python analysis/make_figures.py        # results/figures/*.png
```

`make_figures.py` takes the bar heights from the result files and the 95%
intervals from `results/figures/stats_table.csv`. Pass a folder as argument to
write the figures somewhere else. It needs `numpy` and `matplotlib`.
