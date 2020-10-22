require.config({
    paths: {
        d3: "https://d3js.org/d3.v5.min"
    }
});

require(["d3"], function(d3) {
var globalRegistry = {
    "ids": []
};
globalRegistry["current"] = function() {
    var ids = globalRegistry["ids"];
    for (var i = 0; i < ids.length; ++i) {
        if (globalRegistry[ids[i]].current) {
            return globalRegistry[ids[i]].current.data;
        }
    }
    return null;
}

var allnodes = function(a, nodearray) {
    var elems = a["children"];
    for (var i = 0; i < elems.length; ++i) {
        nodearray.push(elems[i]["name"]);
        allnodes(elems[i], nodearray);
    }
}

var mergenodes = function(a, b) {
    var x = [];
    allnodes(a, x);
    allnodes(b, x);
    var aset = new Set(x);
    var arrayFromSet = [...aset];

    return arrayFromSet;
}


// Dimensions of sunburst.
var width = 680;
var height = 680;
var radius = Math.min(width, height) / 2;

// Breadcrumb dimensions: width, height, spacing, width of tip/tail.
var b = {
    w: 175,
    h: 30,
    s: 3,
    t: 10
};

// Mapping of step names to colors.
var colors = {};

var loadColors = function(names) {

    var colorful = [];
    for (var i = 0; i < 10; i += 2) {
        var j = (i * 172141) % 13371;;
        colorful.push(d3.interpolateSpectral(j / 13371.0));
        colorful.push(d3.schemeCategory10[i]);
        colorful.push(d3.schemeCategory10[i + 1]);
        colorful.push(d3.interpolateViridis(j/13371.0))
    };
    //colorful = colorful.concat(palette);
    var elems = names.slice();
    for (var i = 0; i < elems.length; ++i) {
        var c = d3.color(colorful[i % colorful.length]);
        var y = d3.hsl(c);
        if (y.l > 0.7) {
            y.l = 0.7;
            c = y;
        }
        colors[elems[i]] = c.toString();
    }
}

// Total size of all segments; we set this later, after loading the data.
var partition = d3.partition()
    .size([2 * Math.PI, radius * radius]);

var arc = d3.arc()
    .startAngle(function(d) {
        return d.x0;
    })
    .endAngle(function(d) {
        return d.x1;
    })
    .innerRadius(function(d) {
        return Math.sqrt(d.y0);
    })
    .outerRadius(function(d) {
        return Math.sqrt(d.y1);
    });


// Main function to draw and set up the visualization, once we have the data.
var createVisualization = function(json, mainId, id) {
    var totalSize = 0;

    var idtag = function(klassname) {
        return "#chart" + id + " ." + klassname;
    }

    var idnode = function(nodename) {
        return "#chart" + id + " " + nodename;
    }

    var mainIdNode = function(nodename) {
        return "#main_" + mainId + " " + nodename;
    }

    var vis = d3.select("#chart" + id).append("svg:svg")
        .attr("width", width)
        .attr("height", height)
        .append("svg:g")
        .attr("class", "container")
        .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

    var globalMouseOver = function(d) {
        var ary = globalRegistry["ids"].slice();
        for (var i = 0; i < ary.length; ++i) {
            var entity = globalRegistry[ary[i]];
            if (!(d || d.data)) {
                entity["setCurrent"](null);
            } else {
                entity["setCurrent"](entity["nodesbynames"][d.data.name]);
            }
        }
    };
    var globalMouseOut = function(d) {
        var ary = globalRegistry["ids"].slice();
        for (var i = 0; i < ary.length; ++i) {
            var entity = globalRegistry[ary[i]];
            entity["setCurrent"](null);
        }
    };

    var setCurrentToNull = function(previous) {
        globalRegistry[id].current = null;
        // Hide the breadcrumb trail
        d3.select(mainIdNode(".treelegend")).style("display", "hidden");

        // Deactivate all segments during transition.
        //d3.selectAll(idnode("path")).on("mouseover", null);

        // Transition each segment to full opacity and then reactivate it.
        d3.selectAll(idnode("path"))
            .transition()
            .duration(1000)
            .style("opacity", 1)
        // .on("end", function() {
        //         d3.select(this).on("mouseover", globalMouseOver);
        //       });

        d3.selectAll(idtag("explanation"))
            .style("visibility", "hidden"); // HIDDEN
        updateDescription(globalRegistry.current());
    }

    var copyToClipboard = function(d) {
    	var ary = globalRegistry["ids"].slice();
        for (var i = 0; i < ary.length; ++i) {
            var entity = globalRegistry[ary[i]];
            if (!(d || d.data)) {
                
            } else {
                entity["setCurrent"](entity["nodesbynames"][d.data.name]);

            }
        }
        const svg = d3.select("#chart" + id).select("svg");
    	const str = svg.html();
    	console.log(str);
  		const el = document.createElement('textarea');
  		el.value = str;
  		document.body.appendChild(el);
  		el.select();
  		document.execCommand('copy');
  		document.body.removeChild(el);
	};

    // Fade all but the current sequence, and show it in the breadcrumb trail.
    var setCurrent = function(d) {
        if (d == null) {
            setCurrentToNull(globalRegistry[id].current);
            return;
        }
        globalRegistry[id].current = d;
        var percentage = (100 * d.data.freq).toPrecision(3);
        var percentageString = percentage + "%";
        if (percentage < 0.1) {
            percentageString = "< 0.1%";
        }

        d3.selectAll(idtag("percentage"))
            .text(percentageString);
        d3.selectAll(idtag("category_name")).html("of compounds (<span class=\"catemph\">" + Math.round(d.data.num) + " spectra</span>) in this dataset belong to the category <span class=\"catemph\">" + d.data.name + "</span>.")
        d3.selectAll(idtag("explanation"))
            .style("visibility", "");

        var sequenceArray = d.ancestors().reverse();
        sequenceArray.shift(); // remove root node from the array
        updateBreadcrumbs(sequenceArray, percentageString);

        // Fade all the segments.
        d3.selectAll(idnode("path"))
            .style("opacity", 0.33).transition();

        // Then highlight only those that are an ancestor of the current segment.
        vis.selectAll(idnode("path"))
            .filter(function(node) {
                return (sequenceArray.indexOf(node) >= 0);
            })
            .style("opacity", 1);
        updateDescription(d.data);
    }

    var hideDescription = function() {
        d3.selectAll(idtag("explanation"))
            .style("visibility", "hidden"); // HIDDEN
    }

    var initializeBreadcrumbTrail = function() {
        // Add the svg area.
        var trail = d3.selectAll(idtag("sequence")).append("svg:svg")
            .attr("width", width)
            .attr("height", 50)
            .attr("id", "trail");
        // Add the label at the end, for the percentage.
        trail.append("svg:text")
            .attr("id", "endlabel")
            .style("fill", "#000");
    }

    // Generate a string that describes the points of a breadcrumb polygon.
    var breadcrumbPoints = function(d, i) {
        var points = [];
        points.push("0,0");
        points.push(b.w + ",0");
        points.push(b.w + b.t + "," + (b.h / 2));
        points.push(b.w + "," + b.h);
        points.push("0," + b.h);
        if (i > 0) { // Leftmost breadcrumb; don't include 6th vertex.
            points.push(b.t + "," + (b.h / 2));
        }
        return points.join(" ");
    }

    var updateDescription = function(activeNode) {
        if (activeNode) {
            d3.selectAll(mainIdNode(".description")).html(activeNode.description);
        } else {
            d3.selectAll(mainIdNode(".description")).html("");
        }
    }

    // Update the breadcrumb trail to show the current sequence and percentage.
    var updateBreadcrumbs = function(nodeArray, percentageString) {
        var trail = d3.select(mainIdNode(".treelegend")).selectAll("div.node").data(nodeArray).html(function(x) {
            return x.data.name;
        }).style("background-color", function(x) {
            return colors[x.data.name];
        });
        trail.enter().append("div").attr("class", "node").html(function(x) {
            return x.data.name;
        }).style("background-color", function(x) {
            return colors[x.data.name];
        });
        trail.exit().remove();
        //d3.select("#trail")
        //    .style("visibility", "");
    }

    // Basic setup of page elements.
    initializeBreadcrumbTrail();

    // Bounding circle underneath the sunburst, to make it easier to detect
    // when the mouse leaves the parent g.
    vis.append("svg:circle")
        .attr("r", radius)
        .style("opacity", 0);

    // Turn the data into a d3 hierarchy and calculate the sums.
    var root = d3.hierarchy(json).each(function(x) {
            x.value = x.data.size;
        })
        .sort(function(a, b) {
            return b.value - a.value;
        });

    // For efficiency, filter nodes to keep only those large enough to see.
    var nodes = partition(root).descendants()
        .filter(function(d) {
            return (d.x1 - d.x0 > 0.005); // 0.005 radians = 0.29 degrees
        });

    var nodesbynames = {};
    for (var i = 0; i < nodes.length; ++i) {
        nodesbynames[nodes[i].data.name] = nodes[i];
    }

    globalRegistry[id] = {
        "nodesbynames": nodesbynames,
        "setCurrent": setCurrent,
        "hideDescription": hideDescription,
        "current": null
    };
    globalRegistry["ids"].push(id);

    var path = vis.data([json]).selectAll("path")
        .data(nodes)
        .enter().append("svg:path")
        .attr("display", function(d) {
            return d.depth ? null : "none";
        })
        .attr("d", arc)
        .attr("fill-rule", "evenodd")
        .style("fill", function(d) {
            return colors[d.data.name];
        })
        .style("opacity", 1)
        .on("mouseover", globalMouseOver).on("dblclick",copyToClipboard);

    // Add the mouseleave handler to the bounding circle.
    d3.selectAll(idtag("container")).on("mouseleave", globalMouseOut);

    // Get total size of the tree = value of root node from partition.
    totalSize = path.datum().value;
};
"<CUSTOM-CODE>";
});