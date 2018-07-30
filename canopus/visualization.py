from IPython.display import HTML, display, Javascript
import uuid
from .ontology import Ontology, SiriusWorkspace,CanopusStatistics
import json
import pkg_resources
from pathlib import Path
class CanopusRenderer:
    def __init__(self, workspace, uid=None):
        self.workspace = workspace
        self.uid = uid if uid is not None else str(uuid.uuid4().hex)
        self.treemaps = []
    
    def addTreemap(self,statistics = None):
        if statistics is None:
            statistics = self.workspace.statistics
        self.treemaps.append(statistics)
        
        
    def renderCSS(self):      
        display(HTML(r"""
<style>
#main {
  width: 100%;
}
.sidebar {
  position: absolute;
  margin-right:20px;
  top: 10%;
  bottom:50%;
  right: 0%;
  width:200px;
  height:300px;
  padding: 12px;
  margin-left:1px; margin-right:1px;
  /*background-color: steelblue;*/
  display:inline-block;
  border-radius: 5px;
  border: 2px solid #000000;
  min-width: 100px;
  min-height: 100px;
}
.description {
  font-weight: 600;
  fill: #fff;
  font-size: 14px;
}

.sequence {
  width: 600px;
  height: 70px;
}


.node {
  font-weight: 900;
  color: #fff;
}
div.left {
  float:left;
}
div.right {
  float:right;
}

.chart path, .chart path {
  stroke: #fff;
}
div.chart {
position:relative;
}

h1 {
  text-align: center;
}

.chart > .explanation {
  position: absolute;
  top: 280px;
  left: 250px;
  width: 180px;
  text-align: center;
  z-index: -1;
}


.percentage {
  font-size: 2.5em;
}

.catemph {
  font-weight: bolder;
  font-size:125%;
}

.treelegend {
  min-height:55px;
}
div.node {
  padding: 3px;
  margin-left:1px; margin-right:1px;
  display:inline-block;
  border-radius: 5px;
  border: 2px solid #000000;
}
.sidebar div {
  font-weight: bold;
  fill: #fff;
}
</style>
"""))
        
    def render(self):
        self.renderCSS()
        self.renderHTML()
        self.renderJavascript()
        
    def renderHTML(self):
        charts = "\n".join([
            "<div class=\"chart\" id=\"chart" + self.uid + "_" + str(i) + "\">"
            """
            
                <div class="explanation" style="visibility: hidden;">
                  <span class="percentage"></span><br/>
                  <span class="category_name"></span>
                </div>
              </div>
            """ for (i,_) in enumerate(self.treemaps)])
        display(HTML("<div id=\"main" + "_" +self.uid + "\">" +
"""
      <div class="treelegend"></div>
      <div class="sequence"></div>
""" + charts + """<div class="sidebar">
      <div>
      <p class="description"></p>
      </div>
</div>"""))
        
    def renderJavascript(self):
        customCode = "var json = " + json.dumps(self.workspace.json_treemap()) + """;
var nodes = [];
allnodes(json, nodes);
loadColors(nodes);
        """;
        charts = "\n".join(["createVisualization(" + json.dumps(self.workspace.json_treemap(x)) + ", \"" + self.uid + "\""  + ",\"" + self.uid + "_" + str(i) + "\");" for (i, x) in enumerate(self.treemaps)])
        json_code = pkg_resources.resource_string("canopus.resources", "treemap.js").decode("utf-8")
        display(Javascript(json_code.replace('"<CUSTOM-CODE>";', customCode + "\n" + charts)))
        