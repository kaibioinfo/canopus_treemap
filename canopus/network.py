import xml.etree.ElementTree as ET
import uuid
from IPython.display import HTML, display, Javascript
import pkg_resources
import json
import pandas as pd
import xml.etree.ElementTree as ET
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
    ET.register_namespace('',"http://graphml.graphdrawing.org/xmlns")
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


  def writeWithPieCharts(self, file):
    ns = dict(ml="http://graphml.graphdrawing.org/xmlns")
    xml = ET.parse(self.gnps)
    # add canopus class tag
    c1 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "main class", "attr.type":"string", "for":"node",
     "id":"canopus1"})
    c2 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "upper class", "attr.type":"string", "for":"node",
     "id":"canopus2"})
    c3 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "superclass", "attr.type":"string", "for":"node",
     "id":"canopus3"})
    c4 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "class", "attr.type":"string", "for":"node",
     "id":"canopus4"})
    c5 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "subclass", "attr.type":"string", "for":"node",
     "id":"canopus5"})
    c6 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "siriusQuality", "attr.type":"long", "for":"node",
     "id":"canopus6"})
    c7 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "formula", "attr.type":"string", "for":"node",
     "id":"sirius1"})

    npc_elements = ["canopus_npc0","canopus_npc1","canopus_npc2"]
    npc0 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "npc_pathway", "attr.type":"string", "for":"node",
     "id":"canopus_npc0"})
    npc1 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "npc_superclass", "attr.type":"string", "for":"node",
     "id":"canopus_npc1"})
    npc2 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "npc_class", "attr.type":"string", "for":"node",
     "id":"canopus_npc2"})


    allSecondaryClasses = set([c for c in self.sirius.statistics.secondaryAssignments.values() if c in self.sirius.revmap])
    mappingCl = {("cl%d" % index): name for (index, name) in enumerate(allSecondaryClasses)}
    for key in mappingCl:
        xml.find(".").insert(0, ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": mappingCl[key].name, "attr.type":"float", "for":"node",
     "id":key}))    
    revmap = dict()
    for key in mappingCl:
        revmap[mappingCl[key]] = key
    cother = "clX"
    xml.find(".").insert(0, ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "other", "attr.type":"float", "for":"node",
     "id":"clX"}))    
    xml.find(".").insert(0, c7)
    xml.find(".").insert(0, c6)
    xml.find(".").insert(0, c5)
    xml.find(".").insert(0, c4)
    xml.find(".").insert(0, c3)
    xml.find(".").insert(0, c2)
    xml.find(".").insert(0, c1)
    for c in [npc0,npc1,npc2]:
      xml.find(".").insert(0, c)

    for graph in xml.findall("ml:graph",ns):
      for node in graph.findall("ml:node",ns):
        id = node.attrib["id"]
        if id in self.sirius.compounds:
          c = self.sirius.compounds[id]
          # get quality flag
          flag = 4
          if c.isBadQuality(zodiac=0.8):
            flag = 3
          if c.isBadQuality(peakshape=True):
            flag = 2
          if c.isBadQuality():
            flag = 1
          meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus6"})
          meta.text = str(flag)
          node.append(meta)

          # npc
          npc_list = self.sirius.statistics.npcPrimaryAssignmentVector(c,0.33)
          for i,elem in enumerate(npc_list):
            if elem:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":npc_elements[i]})
              meta.text = elem.name
              node.append(meta)


          if c in self.sirius.statistics.assignments:
            formula = c.formula
            meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"sirius1"})
            meta.text = str(formula)
            node.append(meta)
            klassname = self.sirius.statistics.assignments[c].name
            meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus1"})
            meta.text = klassname
            node.append(meta)
            genus = self.sirius.statistics.assignments[c].classyFireGenus()
            if "superclass" in genus:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus3"})
              meta.text = genus["superclass"].name
              node.append(meta)
            if "class" in genus:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus4"})
              meta.text = genus["class"].name
              node.append(meta)
            if "subclass" in genus:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus5"})
              meta.text = genus["subclass"].name
              node.append(meta)
            klassname2 = self.sirius.statistics.secondaryAssignments[c].name 
            meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus2"})
            meta.text = klassname2
            node.append(meta)
            allSecondaryClasses = set()
            for key in mappingCl:
                index=self.sirius.revmap[mappingCl[key]]
                if c.canopusfp[index]>=0.5:
                    allSecondaryClasses.add(key)

            for key in list(allSecondaryClasses):
                for a in mappingCl[key].ancestors():
                    if a in revmap and revmap[a] in allSecondaryClasses:
                        display("remove %s due to %s" % (a.name, mappingCl[key].name) )
                        allSecondaryClasses.remove(revmap[a])

            weight = 1.0/len(allSecondaryClasses) if allSecondaryClasses else 1.0
            for key in mappingCl:
                meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":key})
                meta.text = str(weight) if (key in allSecondaryClasses) else "0.0"
                node.append(meta)
            otherweight = 0.0 if allSecondaryClasses else 1.0
            meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":cother})
            meta.text = str(otherweight)
            node.append(meta)
    xml.write(file)

  def writeWithPieChartsPreselected(self, file, assignments):
    ns = dict(ml="http://graphml.graphdrawing.org/xmlns")
    xml = ET.parse(self.gnps)
    # add canopus class tag
    c1 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "main class", "attr.type":"string", "for":"node",
     "id":"canopus1"})
    c2 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "upper class", "attr.type":"string", "for":"node",
     "id":"canopus2"})
    c3 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "superclass", "attr.type":"string", "for":"node",
     "id":"canopus3"})
    c4 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "class", "attr.type":"string", "for":"node",
     "id":"canopus4"})
    c5 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "subclass", "attr.type":"string", "for":"node",
     "id":"canopus5"})
    c6 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "siriusQuality", "attr.type":"long", "for":"node",
     "id":"canopus6"})
    c7 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "formula", "attr.type":"string", "for":"node",
     "id":"sirius1"})

    npc_elements = ["canopus_npc0","canopus_npc1","canopus_npc2"]
    npc0 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "npc_pathway", "attr.type":"string", "for":"node",
     "id":"canopus_npc0"})
    npc1 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "npc_superclass", "attr.type":"string", "for":"node",
     "id":"canopus_npc1"})
    npc2 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "npc_class", "attr.type":"string", "for":"node",
     "id":"canopus_npc2"})

    allSecondaryClasses = set()
    for subset in assignments.values():
      for kl in subset:
        allSecondaryClasses.add(kl)
    mappingCl = {("cl%d" % index): name for (index, name) in enumerate(allSecondaryClasses)}
    for key in mappingCl:
        xml.find(".").insert(0, ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": mappingCl[key].name, "attr.type":"float", "for":"node",
     "id":key}))    
    revmap = dict()
    for key in mappingCl:
        revmap[mappingCl[key]] = key
    cother = "clX"
    xml.find(".").insert(0, ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "other", "attr.type":"float", "for":"node",
     "id":"clX"}))    
    xml.find(".").insert(0, c7)
    xml.find(".").insert(0, c6)
    xml.find(".").insert(0, c5)
    xml.find(".").insert(0, c4)
    xml.find(".").insert(0, c3)
    xml.find(".").insert(0, c2)
    xml.find(".").insert(0, c1)
    for c in [npc0,npc1,npc2]:
      xml.find(".").insert(0, c)
    for graph in xml.findall("ml:graph",ns):
      for node in graph.findall("ml:node",ns):
        id = node.attrib["id"]
        if id in self.sirius.compounds:
          c = self.sirius.compounds[id]
          # get quality flag
          flag = 4
          if c.isBadQuality(zodiac=0.8):
            flag = 3
          if c.isBadQuality(peakshape=True):
            flag = 2
          if c.isBadQuality():
            flag = 1
          meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus6"})
          meta.text = str(flag)
          node.append(meta)

          # npc
          npc_list = self.sirius.statistics.npcPrimaryAssignmentVector(c,0.33)
          for i,elem in enumerate(npc_list):
            if elem:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":npc_elements[i]})
              meta.text = elem.name
              node.append(meta)

          if c in self.sirius.statistics.assignments:
            formula = c.formula
            meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"sirius1"})
            meta.text = str(formula)
            node.append(meta)
            klassname = self.sirius.statistics.assignments[c].name
            meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus1"})
            meta.text = klassname
            node.append(meta)
            genus = self.sirius.statistics.assignments[c].classyFireGenus()
            if "superclass" in genus:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus3"})
              meta.text = genus["superclass"].name
              node.append(meta)
            if "class" in genus:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus4"})
              meta.text = genus["class"].name
              node.append(meta)
            if "subclass" in genus:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus5"})
              meta.text = genus["subclass"].name
              node.append(meta)
            klassname2 = self.sirius.statistics.secondaryAssignments[c].name 
            meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus2"})
            meta.text = klassname2
            node.append(meta)
            allSecondaryClasses = assignments[c] if c in assignments else set()

            weight = 1.0/len(allSecondaryClasses) if allSecondaryClasses else 1.0
            for key in mappingCl:
                meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":key})
                meta.text = str(weight) if (key in allSecondaryClasses) else "0.0"
                node.append(meta)
            otherweight = 0.0 if allSecondaryClasses else 1.0
            meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":cother})
            meta.text = str(otherweight)
            node.append(meta)
    xml.write(file)

  def write(self, file, manualAssignment=None):
    ns = dict(ml="http://graphml.graphdrawing.org/xmlns")
    xml = ET.parse(self.gnps)
    # add canopus class tag
    c1 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "main class", "attr.type":"string", "for":"node",
     "id":"canopus1"})
    c2 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "upper class", "attr.type":"string", "for":"node",
     "id":"canopus2"})
    c3 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "superclass", "attr.type":"string", "for":"node",
     "id":"canopus3"})
    c4 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "class", "attr.type":"string", "for":"node",
     "id":"canopus4"})
    c5 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "subclass", "attr.type":"string", "for":"node",
     "id":"canopus5"})
    c6 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "siriusQuality", "attr.type":"long", "for":"node",
     "id":"canopus6"})
    c8 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "classification", "attr.type":"string", "for":"node",
     "id":"canopus7"})
    c7 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "formula", "attr.type":"string", "for":"node",
     "id":"sirius1"})
    npc_elements = ["canopus_npc0","canopus_npc1","canopus_npc2"]
    npc0 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "npc_pathway", "attr.type":"string", "for":"node",
     "id":"canopus_npc0"})
    npc1 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "npc_superclass", "attr.type":"string", "for":"node",
     "id":"canopus_npc1"})
    npc2 = ET.Element("{http://graphml.graphdrawing.org/xmlns}key", {"attr.name": "npc_class", "attr.type":"string", "for":"node",
     "id":"canopus_npc2"})
    xml.find(".").insert(0, c8)
    xml.find(".").insert(0, c7)
    xml.find(".").insert(0, c6)
    xml.find(".").insert(0, c5)
    xml.find(".").insert(0, c4)
    xml.find(".").insert(0, c3)
    xml.find(".").insert(0, c2)
    xml.find(".").insert(0, c1)
    for c in [npc0,npc1,npc2]:
      xml.find(".").insert(0, c)
    for graph in xml.findall("ml:graph",ns):
      for node in graph.findall("ml:node",ns):
        id = node.attrib["id"]
        if id in self.sirius.compounds:
          c = self.sirius.compounds[id]
          # get quality flag
          flag = 4
          if c.isBadQuality(zodiac=0.8):
            flag = 3
          if c.isBadQuality(peakshape=True):
            flag = 2
          if c.isBadQuality():
            flag = 1
          meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus6"})
          meta.text = str(flag)
          node.append(meta)
          formula = c.formula
          meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"sirius1"})
          meta.text = str(formula)
          node.append(meta)

          # npc
          npc_list = self.sirius.statistics.npcPrimaryAssignmentVector(c,0.33)
          for i,elem in enumerate(npc_list):
            if elem:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":npc_elements[i]})
              meta.text = elem.name
              node.append(meta)

          if c in self.sirius.statistics.assignments:
            klassname = self.sirius.statistics.assignments[c].name
            meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus1"})
            meta.text = klassname
            node.append(meta)
            genus = self.sirius.statistics.assignments[c].classyFireGenus()
            if "superclass" in genus:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus3"})
              meta.text = genus["superclass"].name
              node.append(meta)
            if "class" in genus:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus4"})
              meta.text = genus["class"].name
              node.append(meta)
            if "subclass" in genus:
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus5"})
              meta.text = genus["subclass"].name
              node.append(meta)
            klassname2 = self.sirius.statistics.secondaryAssignments[c].name 
            meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus2"})
            meta.text = klassname2
            node.append(meta)
            if manualAssignment is not None and c in manualAssignment:
              klassname2 = manualAssignment[c].name 
              meta = ET.Element("{http://graphml.graphdrawing.org/xmlns}data",{"key":"canopus7"})
              meta.text = klassname2
              node.append(meta)

    xml.write(file)



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
    network.gnps = gnps
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
    self.sirius = sirius
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
  <input type="number" name="networkSel" class="networkSelection" disabled="disabled" min="1" max="100" step="1" value="1"><br />
  <input type="number" name="networkFind" class="networkFind" min="0" max="99999" step="1" value="0">
  <input type="button" name="findNetwork" class="networkFindButton">Find Cluster</button>
  </input>
</form>
      """ + ("<a href=\"\" target=\"_blank\" class=\"clusterLink\">View network in GNPS</a></div>" if self.clusterInfo is not None else "</div>")))

  def renderJavascript(self):
    code = pkg_resources.resource_string("canopus.resources", "network.js").decode("utf-8")
    customCode = self.to_json()
    code = code.replace('"<CUSTOM-CODE>";', customCode + ";\n").replace('<CUSTOM-ID>', self.uid).replace("\"<CUSTOM-ONTOLOGY>\"",self.ontology.to_json())
    display(Javascript(code))

  def __init__(self,uid=None):
    self.nodes = dict()
    self.edges = []
    self.uid = uid if uid is not None else str(uuid.uuid4().hex)