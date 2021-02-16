module GrapetreeCluster

using ArgParse
using CSV
using DataFrames
using LightGraphs
using MetaGraphs
using PhyloNetworks


function arguments(cli)

    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input", "-i"
        help = "Path to Newick tree"
        arg_type = String
        required = true

        "--output", "-o"
        help = "Output path for cluster table"
        arg_type = String
        required = true
    end

    parse_args(cli, s)
end

function main()

    args = arguments(ARGS)

    process(args)

end

function main(cli)

    args = arguments(cli)

    process(args)

end

function process(args)

    graph = load_data(args["input"])

    collapsed_graph = collapse_internal_nodes(graph)

    strains = map(
        node -> get_prop(collapsed_graph, node, :name),
        vertices(collapsed_graph),
    )

    clusters = cluster_graph(collapsed_graph)

    cluster_table = build_cluster_table(strains, clusters)

    write_output(cluster_table, args["output"])

end


function load_data(filepath::String)::MetaGraph

    tree = readTopology(filepath)
    graph = convert_to_lightgraphs(tree)

    graph
end

function convert_to_lightgraphs(tree::HybridNetwork)::MetaGraph

    graph_size = tree.node |> length

    G = MetaGraph(graph_size)

    node_names = Dict()

    for (idx, node) in enumerate(tree.node)
        if isempty(node.name)
            node.name = "hypothetical_$idx"
            is_hypothetical = true
        else
            is_hypothetical = false
        end

        nodename = if isempty(node.name)
            string(idx)
        else
            node.name
        end

        node_names[nodename] = idx

        is_tip = node.leaf
        set_props!(
            G,
            idx,
            Dict(
                :name => nodename,
                :tip => is_tip,
                :hypothetical => is_hypothetical,
            ),
        )
    end

    for edge in tree.edge
        node1, node2 = [n.name for n in edge.node]

        weight = edge.length

        src = node_names[node1]
        dst = node_names[node2]

        add_edge!(G, src, dst, :weight, weight)
    end

    G
end

function collapse_internal_nodes(g::MetaGraph)

    hypotheticals = filter(v -> get_prop(g, v, :hypothetical), vertices(g))
    sort!(hypotheticals, by = x -> length(neighbors(g, x)), rev = true)

    try
        hypothetical = first(hypotheticals)

        trimmed_graph = collapse_internal_node(g, hypothetical)

        return collapse_internal_nodes(trimmed_graph)

    catch BoundsError

        # clean up some self-self edges that occur
        for e in edges(g)
            if e.src == e.dst
                rem_edge!(g, e)
            end
        end

        return g

    end

end

function collapse_internal_node(graph::MetaGraph, vertex::Int)::MetaGraph

    g = deepcopy(graph)
    #g = graph
    # get the neighbours
    neighbours = neighbors(g, vertex)

    # find which neighbour has the shortest distance
    weights =
        [get_prop(g, vertex, neighbour, :weight) for neighbour in neighbours]

    # get the closest, non-self neighbour
    closest_neighbours =
        filter(node -> node != vertex, neighbours[weights.==minimum(weights)])

    closest_neighbour = first(closest_neighbours)

    # replace this hypothetical node with the selected neighbour
    set_props!(g, vertex, props(g, closest_neighbour))

    # transfer edges over
    for neighbour_of_closest in neighbors(g, closest_neighbour)
        w = get_prop(g, closest_neighbour, neighbour_of_closest, :weight)
        add_edge!(g, vertex, neighbour_of_closest, :weight, w)
    end

    # clean up
    rem_vertex!(g, closest_neighbour)

    g
end


function cluster_graph(g::MetaGraph)::Dict{Number,Dict{String,Int}}

    heights = map(e -> get_prop(g, e, :weight), edges(g)) |> sort |> unique

    threshold_clusters = Dict(
        height => cluster_by_delink(g, height) for height in [0, heights...]
    )

    threshold_strain_cluster = Dict(
        height => associate_strain_with_cluster(clusters) for
        (height, clusters) in threshold_clusters
    )

    threshold_strain_cluster
end


function cluster_by_delink(original_graph::MetaGraph, max_distance::Number)

    g = deepcopy(original_graph)

    to_remove = [e for e in edges(g) if get_prop(g, e, :weight) > max_distance]

    for e in to_remove
        rem_edge!(g, e)
    end

    clusters = connected_components(g)

    named_clusters = convert_indices_to_names(g, clusters)

    named_clusters

end


function associate_strain_with_cluster(
    clusters::Array{Array{String,1},1},
)::Dict{String,Int}

    sorted_by_size = sort(clusters, by = length, rev = true)

    strain_cluster = Dict{String,Int}()

    for (cluster_num, members) in enumerate(sorted_by_size)
        for member in members
            strain_cluster[member] = cluster_num

        end
    end

    strain_cluster

end


function convert_indices_to_names(
    graph::MetaGraph,
    clusters::Array{Array{Int,1},1},
)::Array{Array{String,1},1}

    [convert_indices_to_names(graph, cluster) for cluster in clusters]

end


function convert_indices_to_names(
    graph::MetaGraph,
    cluster::Array{Int,1},
)::Array{String,1}

    [get_prop(graph, i, :name) for i in cluster]

end


function build_cluster_table(
    strains::Array{String,1},
    clusters::Dict{Number,Dict{String,Int}},
)::DataFrame

    column_order = clusters |> keys |> collect |> sort

    columns = Dict(
        string(threshold) => cluster_column(threshold, strains, clusters)
        for threshold in keys(clusters)
    )

    df = DataFrame(columns)
    begin
        df[!, "isolate"] = strains
        df
    end

    select(df, ["isolate", map(string, column_order)...])

end

function cluster_column(threshold, strain_order, clusters)::Array{Int,1}

    clusters_at_threshold = clusters[threshold]

    map(x -> clusters_at_threshold[x], strain_order)

end

function write_output(cluster_table::DataFrame, outpath::String)

    CSV.write(outpath, cluster_table; delim = "\t")

end


end # module
