# canopus_treemap
CANOPUS visualization for Jupyter notebook

## usage in Jupyter notebook

```python
import canopus
from canopus.visualization import CanopusRenderer
from canopus.ontology import SiriusWorkspace
# import a workspace with all compounds and annotations
sirius = SiriusWorkspace("/path/to/sirius/workspace")
# display the workspace in Jupyter notebook
r = CanopusRenderer(sirius) 
# visualize all compounds
r.addTreemap()
# visualize just compounds starting with character F
r.addTreemap(sirius.selectByRegexp("^[fF].+"))
# visualize a specific set of compounds
r.addTreemap(sirius.selectByNames(["Flavone", "Nitrendipin", "Cyproterone"]))
# embed visualization in notebook
r.render()
```