"""

	AbqModel(file::AbstractString, inp::File, nodes::Array{Node,1}, minC::Vector, maxC::Vector, dim::Vector,
			 refAxis::AbstractString, defRA::Bool, csys::Array{Int,1}, tol::Float64, defTol::Bool, ecc::Array{String,1},
			 vertices::Dict{AbstractString,Node}, edges::Dict{AbstractString,Array{Node,1}}, faces::Dict{AbstractString,Array{Node,1}},
			 pbcdim::Int, defDim::Bool, eqns::Array{Equation,1}, steps::Array{Step,1})

"""
mutable struct AbqModel
	file::AbstractString
	inp::File
	parts::Dict{String,Part}
	instances::Array{Instance,1}
	nodes::Array{GlobNode,1}
	slaves::Dict{String,Array{Int64,1}}
	minC::Vector
	maxC::Vector
	dim::Vector
	refAxis::AbstractString
	defRA::Bool
	vertexFinder::Bool
	defVF::Bool
	csys::Array{Int,1}
	tol::Float64
	defTol::Bool
	ecc::Array{String,1}
	vertices::Dict{AbstractString,GlobNode}
	edges::Dict{AbstractString,Array{GlobNode,1}}
	faces::Dict{AbstractString,Array{GlobNode,1}}
	pbcdim::Int
	defDim::Bool
	eqns::Array{Equation,1}
	steps::Array{Step,1}
	"""

		AbqModel(file::AbstractString)

	"""
	function AbqModel(file::AbstractString)
		inp = load(file)
		parts, instances, nodes = loadGlobNodes(inp)
		slaves = collect_slaves(inp)
		minC, maxC, dim = getLength(nodes)
		refAxis = "z"
		defRA = true
		vertexFinder = true
		defVF = true
		csys = coords[refAxis]
		tol = 0.001
		defTol = true
		ecc = Array{String,1}()
		vert = Dict{AbstractString,Node}()
		edge = Dict{AbstractString,Array{Node,1}}()
		face = Dict{AbstractString,Array{Node,1}}()
		pbcdim = 3
		defDim = true
		eqns = Array{Equation,1}()
		steps = Array{Step,1}()
		new(file, inp, parts, instances, nodes, slaves, minC, maxC, dim,
			refAxis, defRA, vertexFinder, defVF, csys, tol, defTol, ecc, vert, edge, face, pbcdim, defDim, eqns, steps)
	end
end

"""

	show(io::IO, abq::AbqModel)

"""
function show(io::IO,abq::AbqModel)
	print(io,"AbqModel($(abq.file), $(abq.eqns), $(abq.steps))")
end
 
"""

	coords

Defines constant for easily rotating the boundary conditions with respect to the reference axis.
"""
const coords = Dict("x"=>[1,2,3], "y"=>[2,3,1], "z"=>[3,1,2])

"""

	setRefAxis!(abq::AbqModel, axis::AbstractString)

"""
function setRefAxis!(abq::AbqModel, axis::AbstractString)
	if axis in keys(coords)
		abq.refAxis = axis
		abq.csys = coords[axis]
		abq.defRA = false
		println("Reference axis set to $axis.")
	else
		throw(AxisError)
	end
	return
end

"""

	setPBCdim!(abq::AbqModel, dim::Int)

"""
function setPBCdim!(abq::AbqModel, dim::Int)
	if dim > 0 && dim < 4
		abq.pbcdim = dim
		abq.defDim = false
		println("PBCs are set to $(dim)-dimensional periodicity.")
	else
		throw(DimensionError)
	end
end

"""

	setVertexFinder!(abq::AbqModel, val::Bool)

"""
function setCheckVertices!(abq::AbqModel, val::Bool)
	abq.vertexFinder = val
	abq.defVF = false
	val ? println("Vertices are found automatically") : println("Vertex sets have to be modelled explicitely.")
end

"""

	setTolerance!(abq::AbqModel, newTol::Number)

"""
function setTolerance!(abq::AbqModel, newTol::Number)
	if newTol >= 0
		abq.tol = newTol
		abq.defTol = false
		println("Tolerance set to $newTol.")
	else
		throw(ToleranceError)
	end
	return
end

"""

	setEcceptions!(abq::AbqModel, ecc::Arary{String,1})

"""
function setEcceptions!(abq::AbqModel, ecc::Array{String,1})
	abq.ecc = ecc
	print("Ecceptions set for ")
	n = length(ecc)
	for i = 1:n
		if i == 1
			print(ecc[i])
		elseif i == n
			print(", and $(ecc[i]).\n")
		else
			print(", $(ecc[i])")
		end
	end
	counter = 0
	for i = 1:length(abq.nodes)
		j = i - counter
		if abq.nodes[j].instance in ecc
			deleteat!(abq.nodes,j)
			counter += 1
		end
	end
	abq.minC, abq.maxC, abq.dim = getLength(abq.nodes)
	return
end

"""

	updateNodes!(abq:::AbqModel)

"""
function updateNodes!(abq::AbqModel)
	# Initiate array soon-to-contain the definition of each needed node set for the PBC
	sets = Array{String,1}()
	# 	# Append the nset-definition for each vertex to the array sets
	# 	for v in keys(abq.vertices)
	# 		append!(sets,nset(v,abq.vertices[v].node.num,abq.vertices[v].instance))
	# 	end
	# Append the nset-definition for each edge-node to the array sets
	for e in keys(abq.edges)
		i = 1
		for n in abq.edges[e]
			append!(sets, nset("$(e)-$(i)",n.node.num,n.instance))
			i += 1
		end
	end
	# Append the nset-definition for each edge-node to the array sets
	for f in keys(abq.faces)
		i = 1
		for n in abq.faces[f]
			append!(sets, nset("$(f)-$(i)",n.node.num,n.instance))
			i += 1
		end
	end
	insert!(Line(abq.inp,r"End Assembly"),sets)
	println("Node designation added to input file.")
	return
end	

"""

	LoadCase3D(strain::Matrix{Number}, abq::AbqModel; new=true, name="lc")

Returns a LoadCase defined by the given effective strain tensor strain used for a 3D periodic structure.
"""
function LoadCase3D(strain::Matrix{<:Number}, abq::AbqModel; new=true, name="lc")
	!(length(size(strain))==2 && size(strain, 1)==3 && size(strain,2)==3) && @error "Expected a (3, 3)-array as strain tensor! Got a $(size(strain))-array!"
	boundaries = Array{BoundCon,1}()
	Δu_1 = strain * [1; 0; 0] * abq.dim[abq.csys[1]]
	Δu_2 = strain * [0; 1; 0] * abq.dim[abq.csys[2]]
	Δu_3 = strain * [0; 0; 1] * abq.dim[abq.csys[3]]
	i = 0
	name = new ? "$(name)*" : name
	push!(boundaries, BoundCon("BC-1", new, "SWT", 1, Δu_1[1]))
	push!(boundaries, BoundCon("BC-2", new, "SWT", 2, Δu_1[2]))
	push!(boundaries, BoundCon("BC-3", new, "SWT", 3, Δu_1[3]))
	push!(boundaries, BoundCon("BC-4", new, "NWB", 1, Δu_2[1]))
	push!(boundaries, BoundCon("BC-5", new, "NWB", 2, Δu_2[2]))
	push!(boundaries, BoundCon("BC-6", new, "NWB", 3, Δu_2[3]))
	push!(boundaries, BoundCon("BC-7", new, "SEB", 1, Δu_3[1]))
	push!(boundaries, BoundCon("BC-8", new, "SEB", 2, Δu_3[2]))
	push!(boundaries, BoundCon("BC-9", new, "SEB", 3, Δu_3[3]))
	return LoadCase(name,boundaries)
end

"""

	LoadCase(strain::Matrix{Number}, abq::AbqModel; new=true, name="lc")

Returns a LoadCase defined by the given effective strain tensor strain used for a 3D periodic structure.
"""
function LoadCase(strain::Matrix{<:Number}, abq::AbqModel; new=true, name="lc")
	abq.pbcdim == 1 && (@warn "Function not yet defined for 1D periodicity!"; return LoadCase(name, Array{BoundCon,1}()))
	abq.pbcdim == 2 && (@warn "Function not yet defined for 2D periodicity!"; return LoadCase(name, Array{BoundCon,1}()))
	abq.pbcdim == 3 && return LoadCase3D(strain, abq; new=new, name=name)
end

"""

	LoadCase(name::AbstractString, val::Float64, abq::AbqModel, new::Bool)

"""
function LoadCase(name::AbstractString,val::Float64,abq::AbqModel,new::Bool)
	boundaries = Array{BoundCon,1}()
	bc = loadCases[name]
	i = 0
	name = new ? "$(name)*" : name
	for v in keys(bc)
		dof = bc[v][1]
		disp = bc[v][2]
		for c=1:length(dof)
			i+=1
			if disp[c] == 0
				if new
					push!(boundaries,BoundCon("BC-$(i)",new,v,abq.csys[dof[c]]))
				end
			else
				push!(boundaries,BoundCon("BC-$(i)",new,v,abq.csys[dof[c]],abq.dim[abq.csys[disp[c]]]*val))
			end
		end
	end
	if name == "eps11" || name == "eps11*"
		for n=1:length(abq.faces["B"])
			i+=1
			push!(boundaries,BoundCon("BC-$(i)",new,"B-$(n)",abq.csys[1]))
		end
		for s in ["SB","WB"]
			for n=1:length(abq.edges[s])
				i+=1
				push!(boundaries,BoundCon("BC-$(i)",new,"$(s)-$(n)",abq.csys[1]))
			end
		end
		for n=1:length(abq.faces["T"])
			i+=1
			push!(boundaries,BoundCon("BC-$(i)",new,"T-$(n)",abq.csys[1],abq.dim[abq.csys[1]]*val))
		end
		for s in ["ST","WT"]
			for n=1:length(abq.edges[s])
				i+=1
				push!(boundaries,BoundCon("BC-$(i)",new,"$(s)-$(n)",abq.csys[1],abq.dim[abq.csys[1]]*val))
			end
		end
	elseif name == "eps12" || name == "eps12*"
		for n=1:length(abq.faces["S"])
			i+=1
			push!(boundaries,BoundCon("BC-$(i)",new,"S-$(n)",abq.csys[2],abq.faces["S"][n].coords[abq.csys[1]]*val))
		end
		for s in ["SB","SW","ST","SE"]
			for n=1:length(abq.edges[s])
				i+=1
				push!(boundaries,BoundCon("BC-$(i)",new,"$(s)-$(n)",abq.csys[2],abq.edges[s][n].coords[abq.csys[1]]*val))
			end
		end
	elseif name == "eps13" || name == "eps13*"
		for n=1:length(abq.faces["W"])
			i+=1
			push!(boundaries,BoundCon("BC-$(i)",new,"S-$(n)",abq.csys[3],abq.faces["W"][n].coords[abq.csys[1]]*val))
		end
		for s in ["SW","WT","NW","WB"]
			for n=1:length(abq.edges[s])
				i+=1
				push!(boundaries,BoundCon("BC-$(i)",new,"$(s)-$(n)",abq.csys[3],abq.edges[s][n].coords[abq.csys[1]]*val))
			end
		end
	end
	return LoadCase(name,boundaries)
end

"""

	updatePBC!(abq::AbqModel)

"""
function updatePBC!(abq::AbqModel)
	eqnString = Array{String,1}()
	for e in abq.eqns
		append!(eqnString,generate(e))
	end
	insert!(Line(abq.inp,r"End Assembly"),eqnString)
	println("Periodic boundary conditions added to input file.")
	return
end

"""

	updateSteps!(abq::AbqModel)

"""
function updateSteps!(abq::AbqModel)
	stepString = Array{String,1}()
	steps = deepcopy(abq.steps)
	append!(stepString,["**","** BOUNDARY CONDITIONS","**"])
	iniStep = lbcIni(abq)
	for i in iniStep
		append!(stepString,generate(i))
	end
	for i in steps
		append!(stepString,generate(i))
	end
	append!(Line(abq.inp,length(abq.inp.data)),stepString)
	println("Load boundary conditions added to input file.")
end

"""

	update!(abq::AbqModel)

"""
function update!(abq::AbqModel)
	updateNodes!(abq)
	updatePBC!(abq)
	updateSteps!(abq)
end

"""

	saveInp(abq:::AbqModel, path::AbstractString)

"""
function saveInp(abq::AbqModel,path::AbstractString)
	save(abq.inp,path)
	println("Input file written to $(path)")
	return
end
