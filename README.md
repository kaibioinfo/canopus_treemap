# canopus_treemap
CANOPUS visualization for Jupyter notebook

## install

```
python setup.py install --user
```

## usage in Jupyter notebook

See Tutorial.ipynb for a guide. You can read it from github, but javascript objects won't work. Better is to run it locally:
```
jupyter-notebook Tutorial.ipynb
```

## Changelog

* 27.05.2021 We add support for the [NPC](https://doi.org/10.26434/chemrxiv.12885494.v1) (natural product classifier) prediction. In particular, you can create a summary file with all NPC classes by the following code: 
```python
from canopus import Canopus
C = Canopus(sirius="sirius_projectspace")
C.npcSummary().to_csv("npc_summary.csv")
```

Note that support for NPC predictions exists since several months in CANOPUS. Thus, it is very likely that you do not have to recompute your projectspace. 
