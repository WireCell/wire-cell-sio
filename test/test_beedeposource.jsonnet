local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";

local datadir = std.extVar("datadir");

local depos = g.pnode({
    type: 'BeeDepoSource',
    data: {
        filelist: ["%s/%d/%d-truthDepo.json"%[datadir, n, n] for n in std.range(0,9)],
    },
}, nin=0, nout=1);

local sink = g.pnode({ type: 'DumpDepos' }, nin=1, nout=0);

local graph = g.pipeline([depos, sink]);

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio"],
        apps: ["Pgrapher"]
    },
};

local app = {
    type: "Pgrapher",
    data: {
        edges: g.edges(graph)
    },
};

[cmdline] + g.uses(graph) + [app]
