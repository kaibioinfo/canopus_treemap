import xml.etree.ElementTree as ET
import uuid
from IPython.display import HTML, display, Javascript
import pkg_resources
import json
import pandas as pd
class NetworkEdge:
  def __init__(self, source, target,props):
    self.source = source
    self.target = target
    self.properties = props

  def __getitem__(self,key):
    return self.properties[key]

  def __setitem__(self,key,value):
    self.properties[key] = value

  def reverse(self):
    return NetworkEdge(self.target,self.source,self.properties)

class NetworkNode:
  def __init__(self,nodeId):
    self.edges = []
    self.nodeId = nodeId
    self.gnpsLink = ""
    self.classification = []
    self.singleClass = ""
    self.formula = "?"
    self.componentId = -1

  def connect(self, otherNode):
    p=dict()
    edge = NetworkEdge(self,otherNode,p)
    self.edges.append(edge)
    #otherNode.edges.append(edge.reverse())
    return edge

class MolecularNetwork:

  def parse(gnps):
    ns = dict(ml="http://graphml.graphdrawing.org/xmlns")
    edgeScore = None 
    massDiff = None
    #
    compoundId = None
    mz = None
    gnpsLink = None
    componentId = None
    #
    xml = ET.parse(gnps)
    for key in xml.findall("ml:key",ns):
      n = key.attrib["attr.name"]
      i = key.attrib["id"]
      if n == "EdgeScore":
        edgeScore = i
      elif n == "mass_difference":
        massDiff = i
      elif n == "cluster index":
        compoundId = i
      elif n == "precursor mass":
        mz = i
      elif n == "GNPSLibraryURL":
        gnpsLink = i
      elif n == "componentindex":
        componentId = i

    # now parse all nodes
    network = MolecularNetwork()
    for graph in xml.findall("ml:graph",ns):
      for node in graph.findall("ml:node",ns):
        u = NetworkNode(node.attrib["id"])
        for property in node.findall("ml:data",ns):
          n = property.attrib["key"]
          c = property.text
          if n == compoundId:
            u.compoundId = c
          elif n == mz:
            u.mz = float(c)
          elif n == gnpsLink:
            u.gnpsLink = c
          elif n == componentId:
            u.componentId = c
        network.nodes[u.nodeId] = u
      for edge in graph.findall("ml:edge",ns):
        u = network.nodes[edge.attrib["source"]]
        v = network.nodes[edge.attrib["target"]]
        e = u.connect(v)
        network.edges.append(e)
        for prop in edge.findall("ml:data",ns):
          n = prop.attrib["key"]
          c = prop.text
          if n == edgeScore:
            e["score"] = float(c)
          elif n == massDiff:
            e["delta"] = float(c)
    network.clusterInfo = None
    return network

  def feedClusterInfo(self, file):
    self.clusterInfo = pd.read_csv(file,sep="\t")

  def feedSirius(self, sirius):
    self.ontology = sirius.ontology
    sirius.statistics.assign_most_specific_classes()
    for node in self.nodes.values():
      if node.nodeId in sirius.compounds: 
        compound = sirius.compounds[node.nodeId]
        node.formula = compound.formula
        node.classification = [c.name for c in sirius.statistics.categoriesFor(compound,0.5)]
        node.singleClass = sirius.statistics.assignments[compound].name

  def to_json(self):
    buffer = ["{\"nodes\": ["]
    def prp(k,v=None):
      buffer.append("\"")
      buffer.append(k)
      buffer.append("\": ")
      if v is not None:
        buffer.append(json.dumps(v))
        buffer.append(",")
    def encl(bracket):
      if buffer[len(buffer)-1]==",":
        buffer[len(buffer)-1] = bracket
      else:
        buffer.append(bracket)
    for node in self.nodes.values():
      buffer.append("{")
      prp("nodeId",node.nodeId)
      prp("mz",node.mz)
      prp("url",node.gnpsLink)
      prp("formula",node.formula)
      prp("classification",node.classification)
      prp("category",node.singleClass)
      prp("componentId", node.componentId)
      encl("}")
      buffer.append(",")
    encl("]")
    buffer.append(",\"edges\":[")
    for edge in self.edges:
      buffer.append("{")
      prp("source",edge.source.nodeId)
      prp("target",edge.target.nodeId)
      prp("score",edge["score"])
      prp("delta",edge["delta"])
      encl("}")
      buffer.append(",")
    encl("]")
    if self.clusterInfo is not None:
      buffer.append(",\"networkLinks\": {");  
      for row in self.clusterInfo[self.clusterInfo["componentindex"]!=-1].itertuples():
        prp(row.componentindex, str(row.GNPSLinkout_Network))
      encl("}")
    buffer.append("}")
    return "".join([str(s) for s in buffer])

  def render(self):
    self.renderCSS();
    self.renderHTML();
    self.renderJavascript();

  def renderCSS(self):
    display(HTML(r"""
<style>
svg.network {

}
.glinks line {
  stroke:steelblue;
  stroke-width:1.5px;
}
.gnodes circle {
  fill:steelblue;
  stroke:black;
}
.gnodes circle.identified {
  fill:orange;
}
.gnodes circle.selectedNode, .gnodes circle.identified.selectedNode {
  fill:seagreen;
}
.gnodes circle.highlighted {
  stroke:red;
  stroke-width:2px;
}

form.network {
}
div.networkDescription span.selected {
  text-decoration: underline;
}

div.networkDescription {
  width:400px:
  height:800px;
  border:1px solid;
  float:right;
  overflow: hidden;
}

</style>
"""))

  def renderHTML(self):
    display(HTML("<div id=\"graph" +self.uid + "\"></svg class=\"network\">" + r"""
      <svg class="network"></svg>
      <div class="networkDescription"></div>
      <p>Majority category: <span class='networkMajority'>Unknown</span></p>
<form class="network">
  <input type="radio" name="network" class="networkAll" checked="checked">Display all networks</input>
  <input type="radio" name="network" class="networkSingle">Display subnetwork 
  <input type="number" name="networkSel" class="networkSelection" disabled="disabled" min="1" max="100" step="1" value="1">
  </input>
</form>
      """ + ("<a href=\"\" target=\"_blank\" class=\"clusterLink\">View network in GNPS</a></div>" if self.clusterInfo is not None else "</div>")))

  def renderJavascript(self):
    code = pkg_resources.resource_string("canopus.resources", "network.js").decode("utf-8")
    customCode = self.to_json()
    code = code.replace('"<CUSTOM-CODE>";', customCode + ";\n").replace('<CUSTOM-ID>', self.uid).replace("\"<CUSTOM-ONTOLOGY>\"",self.ontology.to_json())
    with open("fuck.js","w") as f:
      f.write(code)
    display(Javascript(code))

  def __init__(self,uid=None):
    self.nodes = dict()
    self.edges = []
    self.uid = uid if uid is not None else str(uuid.uuid4().hex)