"""

	AbqModel(file::AbstractString, inp::File, nodes::Array{Node,1}, minC::Vector, maxC::Vector, dim::Vector,
			 refAxis::AbstractString, defRA::Bool, csys::Array{Int,1}, tol::Float64, defTol::Bool, exc::Array{String,1},
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
	exc::Array{String,1}
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
		refAxis = "x"
		defRA = true
		vertexFinder = true
		defVF = true
		csys = coords[refAxis]
		tol = 0.001
		defTol = true
		exc = Array{String,1}()
		vert = Dict{AbstractString,Node}()
		edge = Dict{AbstractString,Array{Node,1}}()
		face = Dict{AbstractString,Array{Node,1}}()
		pbcdim = 3
		defDim = true
		eqns = Array{Equation,1}()
		steps = Array{Step,1}()
		new(file, inp, parts, instances, nodes, slaves, minC, maxC, dim,
			refAxis, defRA, vertexFinder, defVF, csys, tol, defTol, exc, vert, edge, face, pbcdim, defDim, eqns, steps)
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
		@info "Reference axis set to $axis."
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
		@info "PBCs are set to $(dim)-dimensional periodicity."
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
	val ? (@info "Vertices are found automatically") : (@info "Vertex sets have to be modelled explicitely.")
end

"""

	setTolerance!(abq::AbqModel, newTol::Number)

"""
function setTolerance!(abq::AbqModel, newTol::Number)
	if newTol >= 0
		abq.tol = newTol
		abq.defTol = false
		@info "Tolerance set to $newTol."
	else
		throw(ToleranceError)
	end
	return
end

"""

	setExceptions!(abq::AbqModel, exc::Arary{String,1})

"""
function setExceptions!(abq::AbqModel, exc::Array{String,1})
	abq.exc = exc
	print("Exceptions set for ")
	n = length(exc)
	for i = 1:n
		if i == 1
			print(exc[i])
		elseif i == n
			print(", and $(exc[i]).\n")
		else
			print(", $(exc[i])")
		end
	end
	counter = 0
	for i = 1:length(abq.nodes)
		j = i - counter
		if abq.nodes[j].instance in exc
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
	if abq.vertexFinder
		for v in keys(abq.vertices)
			append!(sets,nset(v,abq.vertices[v].node.num,abq.vertices[v].instance))
		end
	end
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
	@info "Node designation added to input file."
	return
end

"""

	LoadCase3D(strain::Matrix{Number}, abq::AbqModel;free=true, new=true, name="lc")

Returns a LoadCase defined by the given effective strain tensor strain used for a 3D periodic structure.
"""
function LoadCase3D(strain::Matrix{<:Number}, abq::AbqModel; free=true, new=true, name="lc")
	!(length(size(strain))==2 && size(strain, 1)==3 && size(strain,2)==3) && @error "Expected a (3, 3)-array as strain tensor! Got a $(size(strain))-array!"
	boundaries = Array{BoundCon,1}()
	Δu_1 = strain * [1; 0; 0] * abq.dim[abq.csys[1]]
	Δu_2 = strain * [0; 1; 0] * abq.dim[abq.csys[2]]
	Δu_3 = strain * [0; 0; 1] * abq.dim[abq.csys[3]]
	i = 0
	push!(boundaries, BoundCon("BC-01", new, "SWB", abq.csys[1]))
	push!(boundaries, BoundCon("BC-02", new, "SWB", abq.csys[2]))
	push!(boundaries, BoundCon("BC-03", new, "SWB", abq.csys[3]))
	if free
        Δu_1[1]!=0 && push!(boundaries, BoundCon("BC-04", new, "SWT", abq.csys[1], Δu_1[1]))
        Δu_1[2]!=0 && push!(boundaries, BoundCon("BC-05", new, "SWT", abq.csys[2], Δu_1[2]))
        Δu_1[3]!=0 && push!(boundaries, BoundCon("BC-06", new, "SWT", abq.csys[3], Δu_1[3]))
        Δu_2[1]!=0 && push!(boundaries, BoundCon("BC-07", new, "NWB", abq.csys[1], Δu_2[1]))
        Δu_2[2]!=0 && push!(boundaries, BoundCon("BC-08", new, "NWB", abq.csys[2], Δu_2[2]))
        Δu_2[3]!=0 && push!(boundaries, BoundCon("BC-09", new, "NWB", abq.csys[3], Δu_2[3]))
        Δu_3[1]!=0 && push!(boundaries, BoundCon("BC-10", new, "SEB", abq.csys[1], Δu_3[1]))
        Δu_3[2]!=0 && push!(boundaries, BoundCon("BC-11", new, "SEB", abq.csys[2], Δu_3[2]))
        Δu_3[3]!=0 && push!(boundaries, BoundCon("BC-12", new, "SEB", abq.csys[3], Δu_3[3]))
    else
        push!(boundaries, BoundCon("BC-04", new, "SWT", abq.csys[1], Δu_1[1]))
        push!(boundaries, BoundCon("BC-05", new, "SWT", abq.csys[2], Δu_1[2]))
        push!(boundaries, BoundCon("BC-06", new, "SWT", abq.csys[3], Δu_1[3]))
        push!(boundaries, BoundCon("BC-07", new, "NWB", abq.csys[1], Δu_2[1]))
        push!(boundaries, BoundCon("BC-08", new, "NWB", abq.csys[2], Δu_2[2]))
        push!(boundaries, BoundCon("BC-09", new, "NWB", abq.csys[3], Δu_2[3]))
        push!(boundaries, BoundCon("BC-10", new, "SEB", abq.csys[1], Δu_3[1]))
        push!(boundaries, BoundCon("BC-11", new, "SEB", abq.csys[2], Δu_3[2]))
        push!(boundaries, BoundCon("BC-12", new, "SEB", abq.csys[3], Δu_3[3]))
    end
	return LoadCase(name,boundaries)
end

"""

	LoadCase(strain::Matrix{Number}, abq::AbqModel; new=true, name="lc")

Returns a LoadCase defined by the given effective strain tensor strain used for a 3D periodic structure.
"""
function LoadCase(strain::Matrix{<:Number}, abq::AbqModel;free=free, new=true, name="lc")
	abq.pbcdim == 1 && (@warn "Function not yet defined for 1D periodicity!"; return LoadCase(name, Array{BoundCon,1}()))
	abq.pbcdim == 2 && (@warn "Function not yet defined for 2D periodicity!"; return LoadCase(name, Array{BoundCon,1}()))
	abq.pbcdim == 3 && return LoadCase3D(strain, abq;free=free, new=new, name=name)
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
	@info "Periodic boundary conditions added to input file."
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
	@info "Load boundary conditions added to input file."
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
	@info "Input file written to $(path)"
	return
end
