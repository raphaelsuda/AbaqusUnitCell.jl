#=
	Constant for defining location of vertices
	Each array specifies the positioning of the
	in terms of relative I-II-III-coordinates.
=#
const vertices = Dict("SWB"=>[0,0,0],
					  "SWT"=>[1,0,0],
					  "NWB"=>[0,1,0],
					  "NWT"=>[1,1,0],
					  "SEB"=>[0,0,1],
					  "SET"=>[1,0,1],
					  "NEB"=>[0,1,1],
					  "NET"=>[1,1,1])

#=
	Constant for defining location of edges
	Each array specifies the positioning of the
	in terms of relative I-II-III-coordinates.
	The first number indicates the axis in which
	the edge is oriented, the remaining numbers
	define the coordinates of the edge normal to
	the orientation axis.
=#
const edges = Dict("SW"=>[1,0,0],
				   "SE"=>[1,0,1],
				   "NW"=>[1,1,0],
				   "NE"=>[1,1,1],
				   "WB"=>[2,0,0],
				   "WT"=>[2,1,0],
				   "EB"=>[2,0,1],
				   "ET"=>[2,1,1],
				   "SB"=>[3,0,0],
				   "ST"=>[3,1,0],
				   "NB"=>[3,0,1],
				   "NT"=>[3,1,1])
#=
	Define which coordinates are defined by
	the second and third number of the arrays
	in edges.
=#
const edgeCS = [2 3;1 3;1 2]

#=
	Constant for defining location of faces
	Each array specifies the positioning of the
	in terms of relative I-II-III-coordinates.
	The first number indicates the normal axis of
	the face, while the second number defines the
	faces coordinate in that direction.
=#
const faces = Dict("B"=>[1,0],
				   "T"=>[1,1],
				   "S"=>[2,0],
				   "N"=>[2,1],
				   "W"=>[3,0],
				   "E"=>[3,1])

const faceCS = [2 3;1 3;1 2]

"""

	isEqual(num1::Number,num2::Number,tolerance::Number)

"""
function isEqual(num1::Number,num2::Number,tolerance::Number)
	val = num1 > num2 - tolerance &&
		  num1 < num2 + tolerance
	return val
end

"""

	checkNum(d::Dict,k1::AbstractString,k2::AbstractString)

"""
function checkNum(d::Dict,k1::AbstractString, k2::AbstractString)
	if length(d[k1]) != length(d[k2])
		println("Node numbers of $k1 and $k2 are not equal!")
	end
	return
end

"""

	sortNodes(nodes::Array{GlobNode,1})

Function for sorting of an array containing nodes, according to their coordinates.
Sorting considers x-coordinates first, y-coordinates second, and z-coordinates last.
"""
function sortNodes(nodes::Array{GlobNode,1})
	# Array with weighting factors for sorting of nodes
	weights = [1e11 1e5 1]
	# Initiate array, soon-to-hold the values decisive for sorting
	sortVals = Array{Float64,1}()
	# Fill array sortVals with values, calculated by coordinates
	for n in nodes
		sV = weights * n.node.coords
		append!(sortVals,sV)
	end
	# Obtain array containing the order of the sorted elements
	sortOrder = sortperm(sortVals)
	# Return sorted array of nodes
	return nodes[sortOrder]
end


"""

	findVertex(abq::AbqModel, name::AbstractString)

"""
function findVertex(abq::AbqModel, name::AbstractString)
	nodeBool = map(abq.nodes) do n
		!(n.instance in abq.ecc) && # Is instance containing node in ecceptions?
		isEqual(n.node.coords[abq.csys[1]], abq.minC[abq.csys[1]] + vertices[name][1] * abq.dim[abq.csys[1]], abq.tol) &&
		isEqual(n.node.coords[abq.csys[2]], abq.minC[abq.csys[2]] + vertices[name][2] * abq.dim[abq.csys[2]], abq.tol) &&
		isEqual(n.node.coords[abq.csys[3]], abq.minC[abq.csys[3]] + vertices[name][3] * abq.dim[abq.csys[3]], abq.tol)
	end
	if length(abq.nodes[nodeBool]) != 1
		throw(ToleranceError)
	end
	return abq.nodes[nodeBool][1]
end

"""

	findVertices!(abq::AbqModel)

"""
function findVertices!(abq::AbqModel)
	for v in keys(vertices)
		abq.vertices[v] = findVertex(abq, v)
	end
	println("Vertices written to AbqModel.")
	return
end

"""

	findEdge(abq::AbqModel, name::AbstractString)

"""
function findEdge(abq::AbqModel, name::AbstractString)
	axis1 = abq.csys[edges[name][1]]
	axis2 = abq.csys[edgeCS[edges[name][1],1]]
	axis3 = abq.csys[edgeCS[edges[name][1],2]]
	nodeBool = map(abq.nodes) do n
		!(n.instance in abq.ecc) && # Is instance containing node in ecceptions?
		n.node.coords[axis1] > abq.minC[axis1] + abq.tol &&
		n.node.coords[axis1] < abq.minC[axis1] + abq.dim[axis1] - abq.tol &&
		isEqual(n.node.coords[axis2], abq.minC[axis2] + edges[name][2]*abq.dim[axis2], abq.tol) &&
		isEqual(n.node.coords[axis3], abq.minC[axis3] + edges[name][3]*abq.dim[axis3], abq.tol)
	end
	return abq.nodes[nodeBool]
end

"""

	findEdges!(abq::AbqModel)

"""
function findEdges!(abq::AbqModel)
	for e in keys(edges)
		abq.edges[e] = sortNodes(findEdge(abq, e))
	end
	# Check if corresponding edges have equal node numbers
	checkNum(abq.edges,"SW","NW")
	checkNum(abq.edges,"SW","SE")
	checkNum(abq.edges,"SW","NE")
	checkNum(abq.edges,"WB","EB")
	checkNum(abq.edges,"WT","ET")
	checkNum(abq.edges,"SB","NB")
	checkNum(abq.edges,"ST","NT")
	println("Edges written to AbqModel.")
	return
end


"""

	findFace(abq::AbqModel, name::AbstractString)

"""
function findFace(abq::AbqModel, name::AbstractString)
	axis1 = abq.csys[faces[name][1]]
	axis2 = abq.csys[faceCS[faces[name][1],1]]
	axis3 = abq.csys[faceCS[faces[name][1],2]]
	nodeBool = map(abq.nodes) do n
		!(n.instance in abq.ecc) && # Is instance containing node in ecceptions?
		n.node.coords[axis2] > abq.minC[axis2] + abq.tol &&
		n.node.coords[axis2] < abq.minC[axis2] + abq.dim[axis2] - abq.tol &&
		n.node.coords[axis3] > abq.minC[axis3] + abq.tol &&
		n.node.coords[axis3] < abq.minC[axis3] + abq.dim[axis3] - abq.tol &&
		isEqual(n.node.coords[axis1], abq.minC[axis1] + faces[name][2]*abq.dim[axis1], abq.tol)
	end
	return abq.nodes[nodeBool]
end

"""

	findFaces!(abq::AbqModel)

"""
function findFaces!(abq::AbqModel)
	for f in keys(faces)
		abq.faces[f] = sortNodes(findFace(abq, f))
	end
	# Check if faces West and East have equal number of nodes
	checkNum(abq.faces,"W","E")
	# Check if faces South and North have equal number of nodes
	checkNum(abq.faces,"S","N")
	println("Faces written to AbqModel.")
	return
end

"""

	nset(name::AbstractString, node::Int, inst::AbstractString)

Generate an array containing the definition of a node set in a Abaqus input file.
"""
function nset(name::AbstractString, node::Int, inst::AbstractString)
	return ["*Nset, nset=$(name), instance=$(inst)",
			" $(node),"]
end


"""

	nodeDesignation!(abq::AbqModel)

"""
function nodeDesignation!(abq::AbqModel)
	abq.defRA && @info "Reference axis not explicitly set. Default $(abq.refAxis) is used."
	abq.defVF && @info "Vertex finder not explicitely set. Vertices are being found automatically by default."
	abq.vertexFinder && findVertices!(abq)
	findEdges!(abq)
	findFaces!(abq)
	println("Node designation written to AbqModel.")
	return
end


