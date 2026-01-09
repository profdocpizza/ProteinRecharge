# ProteinRecharge

<p align="center">
  <img src="images/ProteinRecharge_transparent.png" alt="ProteinRecharge logo" width="500"/>
</p>

ProteinRecharge is a small utility for recharging protein surfaces after binder design pipelines produce inappropriately charged binders, using LigandMPNN outputs.

## Quick start ✅

1. Clone this repository and install LigandMPNN:

```bash
git clone https://github.com/profdocpizza/ProteinRecharge.git
cd ProteinRecharge

```

2. Follow [LigandMPNN](https://github.com/dauparas/LigandMPNN/) instructions to install it and download weights 

3. Edit `examples/binder_recharge.yaml` and set the following fields for your environment:

```yaml
ligandmpnn_run_py: "/path/to/LigandMPNN/run.py"   # path to the LigandMPNN runner
ligandmpnn_conda_env: "ligandmpnn_env"          # or "" to use the current env
```

4. Run the pipeline:

```bash
python recharge.py examples/binder_recharge.yaml
```

## Contributing

PRs welcome — please include a short description and a reproducible example when opening issues.
