# ProteinRecharge

ProteinRecharge is a small utility for recharging protein surfaces after binder design pipelines produce inappropriately charged binders, using LigandMPNN outputs.

## Quick start ✅

1. Clone this repository and install LigandMPNN:

```bash
git clone <repo-url>
cd ProteinRecharge   # or `cd recharge` if you haven't renamed the repo folder locally
# Follow LigandMPNN instructions and install it, e.g.:
# https://github.com/dauparas/LigandMPNN/
```

2. Download pretrained LigandMPNN model parameters into the `model_params/` folder of your LigandMPNN installation (see the LigandMPNN README).

3. Edit `examples/binder_recharge.yaml` and set the following fields for your environment:

```yaml
ligandmpnn_run_py: "/path/to/LigandMPNN/run.py"   # path to the LigandMPNN runner
ligandmpnn_conda_env: "ligandmpnn_env"          # or "" to use the current env
```

4. Run the pipeline:

```bash
python recharge.py examples/binder_recharge.yaml
```

Notes
- The script `recharge.py` invokes LigandMPNN using the directory that contains `run.py` as the working directory, so relative model paths (e.g., `./model_params/...`) will resolve correctly.
- If you prefer, you can run `LigandMPNN` in an activated conda env and set `ligandmpnn_conda_env` to `""` in the config.

## Outputs & logs

- Results are written to the `out_folder` specified in the config. The `outputs/` folder is git-ignored and should not be committed.

## Contributing

PRs welcome — please include a short description and a reproducible example when opening issues.
