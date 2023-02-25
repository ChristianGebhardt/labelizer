from pathlib import Path

from labelizer import Fluorophore, Labelizer, label_score_model_paper

cwd = Path(__file__).parent


def main():
    workdir = Path(__file__).parent.joinpath("workdir")
    workdir.mkdir(exist_ok=True)
    labelizer = Labelizer(
        "s",
        str(cwd.joinpath("1DDB/1DDB-39.pdb")),
        None,
        ["A"],
        [],
        label_score_dict=label_score_model_paper,
        workdir=str(workdir),
        prot1_cs=[str(cwd.joinpath("1DDB/1DDB-39_cs.pdb"))],
        prot2_cs=None,
        save_pdb=True,
        cons_serv_add="",
    )
    fluorophore_1 = Fluorophore()
    fluorophore_1.load_from_drive("Cy3")
    labelizer.set_fluorophore(fluorophore_1, None)
    print("Calculating parameter score")
    labelizer.calc_parameter_score()
    print("Calculating labeling score")
    labelizer.calc_labeling_score()
    print("Calculating done")


if __name__ == "__main__":
    main()
