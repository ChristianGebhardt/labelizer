# Labelizer

See our Publication [Labelizer: systematic selection of label sites on proteins (pending)](#)

## Usage

```python
from pathlib import Path

from labelizer import Fluorophore, Labelizer, label_score_model_paper

workdir = Path(__file__).parent.joinpath("workdir")
workdir.mkdir(exist_ok=True)
labelizer = Labelizer(
    "s",
    str(Path("examples/1DDB/1DDB-39.pdb")),
    None,
    ["A"],
    [],
    label_score_dict=label_score_model_paper,
    workdir=str(workdir),
    prot1_cs=[str(Path("examples/1DDB/1DDB-39_cs.pdb"))],
    prot2_cs=None,
    save_pdb=True,
    cons_serv_add="",
)
fluorophore_1 = Fluorophore()
fluorophore_1.load_from_drive("Cy3")
labelizer.set_fluorophore(fluorophore_1, None)
print(f"Calculating parameter score")
labelizer.calc_parameter_score()
print(f"Calculating labeling score")
labelizer.calc_labeling_score()
print(f"Calculating done")
```