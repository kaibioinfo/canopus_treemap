"use strict";
require.config({
    paths: {
        d3: "https://d3js.org/d3.v5.min"
    }
});

require(["d3"], function(d3) {

    function findSubnetworks(graph) {
        let subnetworks = [];
        const node2subnetwork = {};
        const name2nodes = {};

        for(var i=0; i < graph.nodes.length; ++i) {
            graph.nodes[i].neighbours = [];
            graph.nodes[i].classificationSet = d3.set(graph.nodes[i].classification);
            name2nodes[graph.nodes[i].nodeId] = graph.nodes[i];
        }

        for(var i=0; i < graph.edges.length; ++i) {
            const edge = graph.edges[i];
            name2nodes[edge.source].neighbours.push(name2nodes[edge.target]);
            name2nodes[edge.target].neighbours.push(name2nodes[edge.source]);
        }

        function visit(node,subnetwork) {
            subnetwork.nodes.push(node);
            node.subnetwork = subnetwork;
            for (var j=0; j < node.neighbours.length; ++j) {
                const n = node.neighbours[j];
                const c = n.subnetwork;
                n.subnetwork = subnetwork;
                if (c != subnetwork) {
                    visit(n,subnetwork);
                }
            }
        }
        function newNetwork() {
            const network = {};
            network.nodes = [];
            network.edges = [];
            subnetworks.push(network);
            return network;
        }
        var network;
        for(var i=0; i < graph.nodes.length; ++i) {
            if (!graph.nodes[i].subnetwork) {
                network = newNetwork();
                visit(graph.nodes[i],network);
            }
        }
        subnetworks = subnetworks.filter((e,i)=>e.nodes.length>1);
        for(var i=0; i < graph.edges.length; ++i) {
            const e = graph.edges[i];
            name2nodes[e.source].subnetwork.edges.push(e);
        }
        subnetworks.sort((u,v)=>v.nodes.length-u.nodes.length);
        return subnetworks;
    }

    const id = "graph<CUSTOM-ID>";
    const graph = "<CUSTOM-CODE>";

    const ontology = "<CUSTOM-ONTOLOGY>";

    const rootDiv = d3.select("div#"+id);

    const majorityNet = rootDiv.select(".networkMajority")

    const clusterLink = rootDiv.select("a.clusterLink")
    clusterLink.style("visibility","hidden")

    const ontologyTree = ontology["Chemical entities"];
    for (var key in ontology) {
        if (ontology.hasOwnProperty(key)) {
            ontology[key].children = [];
            ontology[key].name = key;
            ontology[key].level = 0;
        }
    }
    var ontologyRoot;
    for (var key in ontology) {
        if (ontology.hasOwnProperty(key)) {
            if (ontology[key].parent.length>0) {
                ontology[key].parentNode = ontology[ontology[key].parent];
                ontology[key].parentNode.children.push(ontology[key]);
            } else ontologyRoot = ontology[key];
        }
    }
    var maxOntologyLevel = 1;
    // add level annotation
    (function(){
        var node = ontologyRoot;
        node.level = 1;
        const stack = [];
        for (var i=0; i < node.children.length; ++i) {
            stack.push(node.children[i]);
        }
        while (stack.length>0) {
            const u = stack.pop();
            u.level = u.parentNode.level + 1;
            maxOntologyLevel = Math.max(u.level,maxOntologyLevel);
            for (var i=0; i < u.children.length; ++i) {
                stack.push(u.children[i]);
            }
        }
    })();

    var selectedClass = null;
    function selectKlass(name) {
        if(selectedClass){
            textbox.select("span[data-id=\"" + selectedClass + "\"]").classed("selected",false);
        }
        selectedClass = name;
        textbox.select("span[data-id=\"" + selectedClass + "\"]").classed("selected",true);
        var nodes = rootDiv.select("svg").selectAll("circle");
        nodes.classed("highlighted",false);
        for (var i=0; i < currentNetwork.nodes.length; ++i) {
            if (currentNetwork.nodes[i].classificationSet.has(selectedClass)) {
                nodes.filter("circle[data-id=\"" + currentNetwork.nodes[i].nodeId + "\"]").classed("highlighted",true);
            }
        }
        majorityNet.classed("selected", name == majorityNet.text());

    }
    function writeList(categories, ht) {
        ht=ht.append("ul")
        function p(node, html) {
            for (var i=0; i < node.children.length; ++i) {
                const c = node.children[i];
                if( categories.has(c.name)) {
                    const t=html.append("li");
                    t.append("span").attr("data-id",c.name).text(c.name).on("click",function(){
                        selectKlass(c.name);
                    });
                    p(c,t.append("ul"));
                }
            }
        }
        p(ontologyTree,ht);
    }

    const allNetworks = findSubnetworks(graph);
    const textbox = rootDiv.select("div.networkDescription");
    const w = rootDiv.node().parentElement.clientWidth;
    const h = w/2.0;
    textbox.attr("width",400);

    var selectedNode = null;

    function displayNode(d) {
        textbox.html("");
        if (!d) return;
        const h = textbox.append("h3")
        if (d.url.length>0) {
            h.append("a").attr("href",d.url).text(d.nodeId + " " + d.formula);
        } else {
            h.text(d.nodeId + " " + d.formula);
        }

        textbox.append("br");
        writeList(d.classificationSet, textbox);
    }

    var currentNetwork = {};

    const networkVisualizer = {
        create: function(Nodes, Edges,connected){
            currentNetwork.nodes = Nodes;
            currentNetwork.edges = Edges;
            const W = w-450;
            this.svg = rootDiv.select("svg");
            this.svg.attr("width",W).attr("height",h);
            const r = 3;
            this.sim = d3.forceSimulation(Nodes)
                .force("link", d3.forceLink(Edges).id(d=>d.nodeId))
                .force("charge", d3.forceManyBody().strength(connected ? -20 : -40))
                .force("bounding", d3.forceCollide(3*r))
                .on("tick",ticked);
            if (connected) {
                //this.sim.force("center",d3.forceCenter(W/2.0,h/2.0));
                this.sim.force("x",d3.forceX(W/2.0));
                this.sim.force("y",d3.forceY(h/2.0));
            } else {
                this.sim.force("x",d3.forceX(W/2.0));
                this.sim.force("y",d3.forceY(h/2.0));
            }

            const container = this.svg.append("g");
            const edgeContainer = container.append("g").attr("class","glinks");
            const nodeContainer = container.append("g").attr("class","gnodes");
            const links = edgeContainer.selectAll("line").data(Edges).enter().append("line");
            const nodes = nodeContainer.selectAll("circle").data(Nodes).enter().append("circle").attr("r",r).classed("identified",d=>d.url.length>0).attr("data-id",d=>d.nodeId)
                .on("mouseover", d=>displayNode(d)).on("mouseout",d=>displayNode(selectedNode)).on("click",function(d){
                    if (selectedNode) {
                        nodeContainer.select("circle[data-id=\""+selectedNode.nodeId+"\"]").classed("selectedNode",false).attr("r",r);
                    }
                    selectedNode = d;
                    displayNode(selectedNode);
                    nodeContainer.select("circle[data-id=\""+selectedNode.nodeId+"\"]").classed("selectedNode",true).attr("r",r*2);
                });

            function ticked() {
                links.attr("x1",d=>d.source.x);
                links.attr("x2",d=>d.target.x);
                links.attr("y1",d=>d.source.y);
                links.attr("y2",d=>d.target.y);
                nodes.attr("cx",d=>d.x).attr("cy",d=>d.y);
            }

            // find majority class
            const allClasses = d3.map();
            for (var c=0; c < currentNetwork.nodes.length; ++c) {
                const node = currentNetwork.nodes[c];
                for (var d = 0; d < node.classification.length; ++d) {
                    if (!allClasses.get(node.classification[d])) {
                        allClasses.set(node.classification[d], 1);
                    } else {
                        allClasses.set(node.classification[d], allClasses.get(node.classification[d]) + 1);
                    }
                }
            }
            var majority = ["unknown",0,0.0];
            var highest = ["unknown",0];
            allClasses.each((value,key)=>{
                const weighted = value + (ontology[key].level/maxOntologyLevel);
                if (weighted > highest[1]) {
                    highest[1] = weighted;
                    highest[0] = key;
                }
                if (value >= currentNetwork.nodes.length/2.0 && (ontology[key].level > majority[1] || (ontology[key].level == majority[1] && weighted>majority[2]))) {
                    majority[0] = key;
                    majority[1] = ontology[key].level; 
                    majority[2] = weighted;
                }
            });
            if (majority[1]>0) {
                majorityNet.text(majority[0]);
                selectKlass(majority[0]);
            } else {
                majorityNet.text(highest[0]);
                selectKlass(majority[0]);
            }

            if (connected) {
                for (var n=0; n < currentNetwork.nodes.length; ++n) {
                    const cid = currentNetwork.nodes[n].componentId;
                    if (cid != "-1") {
                        clusterLink.style("visibility","hidden");
                        var link = graph.networkLinks;
                        if (link) {
                            link = link[cid]
                            if (link) {
                                clusterLink.style("visibility","visible");
                                clusterLink.attr("href", link)
                            }
                        }
                    }
                }
            }
        },
        update: function(Nodes,Edges,connected) {
            if (this.sim) {
                this.exit();
            }
            this.create(Nodes,Edges,connected);
        },
        exit: function() {
            this.sim.on("tick",null);
            this.sim.stop();
            this.svg.html("");
            this.sim = null;
        }
    };

    var selectedNetworkID = 0;
    var isDisplayingSubnetwork = false;

    majorityNet.on("click", function(){
        selectKlass(majorityNet.text());
    });


    const form = rootDiv.select("form");
    form.select("input.networkAll").on("change",function(){
        isDisplayingSubnetwork = !this.checked;
        updateRendering();
    });
    const slider = form.select("input.networkSelection");
    slider.attr("max",allNetworks.length).on("change",function(e){
        selectedNetworkID = parseInt(this.value)-1;
        updateRendering();
    });
    form.select("input.networkSingle").on("change",function(){
        isDisplayingSubnetwork = this.checked;
        updateRendering();
    })

    form.select("input.networkFindButton").on("click", function() {
        var clusterid = parseInt(d3.select("input.networkFind").property("value"));
        for (var j=0; j < allNetworks.length; ++j) {
            for (var k=0; k < allNetworks[j].nodes.length; ++k) {
                if (allNetworks[j].nodes[k].nodeId == clusterid) {
                    slider.property("value", j);
                    selectedNetworkID = j;
                    isDisplayingSubnetwork = true;
                    updateRendering();
                    return;
                }
            }
        }
    });

    function updateRendering() {
        if (isDisplayingSubnetwork) {
            slider.attr("disabled",null);
            networkVisualizer.update(allNetworks[selectedNetworkID].nodes, allNetworks[selectedNetworkID].edges,true );
        } else {
            slider.attr("disabled",true);
            networkVisualizer.update(graph.nodes, graph.edges,false);
        }
    }

    networkVisualizer.create(graph.nodes,graph.edges,false);
    isDisplayingSubnetwork = true;
    updateRendering();
});