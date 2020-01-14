module AbaqusUnitCell

export AbqModel,
	   pbc!,
	   nodeDesignation!,
	   setRefAxis!,
	   setTolerance!,
	   setEcceptions!,
	   setPBCdim!,
	   update!,
	   saveInp,
	   LoadCase,
	   Output,
	   addStep!,
	   +

using StringArrayEditor
import StringArrayEditor: File
import Base:show
import Base:+

"""

	Node(num::Int, coords::Vector)

"""
mutable struct Node
	num::Int
	coords::Vector
end

"""

	GlobNode(glob::Int,instance::String, node::Node)

"""
mutable struct GlobNode
	glob::Int
	instance::String
	node::Node
end

"""

Instance(part::AbstractString, name::AbstractString, trans::Vector, rot_start::Vector, rot_end::Vector, angle_deg::Number, node::Array{Nodes,1})

"""
mutable struct Instance
	part::AbstractString
	name::AbstractString
	trans::Vector
	rot_start::Vector
	rot_end::Vector
	angle_deg::Number
	nodes::Array{Node,1}
	"""

		Instance(part::AbstractString, name::AbstractString)

	"""
	function Instance(part::AbstractString,name::AbstractString)
		new(part, name, [0,0,0],[0,0,0],[1,0,0],0,Node[])
	end
	"""

		Instance(part::AbstractString, name::AbstractString, trans::Vector)

	"""
	function Instance(part::AbstractString,name::AbstractString,trans::Vector)
		new(part, name, trans, [0,0,0], [1,0,0], 0,Node[])
	end
	"""

		Instance(part::AbstractString, name::AbstractString, trans::Vector, rot_start::Vector, rot_end::Vector, angle_deg::Number)

	"""
	function Instance(part::AbstractString,name::AbstractString,trans::Vector, rot_start::Vector, rot_end::Vector, angle_deg::Number)
		new(part, name, trans, rot_start, rot_end, angle_deg, Node[])
	end
end

"""

	Equation(num::Int, nodes::Array{String,1}, dof::Array{Int,1}, coef::Array{Float64,1})

"""
mutable struct Equation
	num::Int
	nodes::Array{String,1}
	dof::Array{Int,1}
	coef::Array{Float64,1}
end

"""

	BoundCon(name::AbstractString, type::AbstractString, new::Bool, node::AbstractString, dof::Int, disp::Float64)

"""
mutable struct BoundCon
	name::AbstractString
	type::AbstractString
	new::Bool
	node::AbstractString
	dof::Int
	disp::Float64
	"""

		BoundCon(name::AbstractString, node::AbstractString, dof::Int)

	"""
	function BoundCon(name::AbstractString,node::AbstractString,dof::Int)
		new(name, "Displacement/Rotation", false, node, dof, 0.0)
	end
	"""

		BoundCon(name::AbstractString, node::AbstractString, dof::Int, disp::Float64)

	"""
	function BoundCon(name::AbstractString,node::AbstractString,dof::Int,disp::Float64)
		new(name, "Displacement/Rotation", false, node, dof, disp)
	end
	"""

		BoundCon(name::AbstractString, new::Bool, node::AbstractString, dof::Int)

	"""
	function BoundCon(name::AbstractString,new::Bool,node::AbstractString,dof::Int)
		new(name, "Displacement/Rotation", new, node, dof, 0.0)
	end
	"""

		BoundCon(name::AbstractString, new::Bool, node::AbstractString, dof::Int, disp::Float64)

	"""
	function BoundCon(name::AbstractString,new::Bool,node::AbstractString,dof::Int,disp::Float64)
		new(name, "Displacement/Rotation", new, node, dof, disp)
	end
end

"""

	+(a::BoundCon, b::BoundCon)

"""
function +(a::BoundCon,b::BoundCon)
	if a.node == b.node && a.dof == b.dof
		BoundCon(a.name,a.new,a.node,a.dof,a.disp+b.disp)
	else
		throw(NodeError)
	end
end

"""

	+(a::Array{BoundCon,1}, b::Array{BoundCon,1})

"""
function +(a::Array{BoundCon,1},b::Array{BoundCon,1})
	bc = deepcopy(a)
	for bcb in b
		bcBool = map(bc) do bca
			bca.node == bcb.node && bca.dof == bcb.dof
		end
		if isempty(bc[bcBool])
			push!(bc,bcb)
		else
			bc[bcBool] = [bc[bcBool][1] + bcb]
		end
	end
	return bc
end

"""

	LoadCase(name::AbstractString, bound::Array{BoundCon,1})

"""
mutable struct LoadCase
	name::AbstractString
	bound::Array{BoundCon,1}
end

"""

	+(a::LoadCase, b::LoadCase)

"""
function +(a::LoadCase,b::LoadCase)
	LoadCase(string(a.name,"+",b.name),a.bound+b.bound)
end

"""

	Output(type::AbstractString, sign::AbstractString)

"""
mutable struct Output
	type::AbstractString
	sign::AbstractString
end

"""

	Step(name::AbstractString, nl::Bool, inc::Int, stab::Float64, allsdtol::Float64, iStart::Float64, iTot::Float64, iMin::Float64, iMax::Float64, nImax::Int, loadCase::LoadCase, output::Array{Output,1})

"""
mutable struct Step
	name::AbstractString
	nl::Bool
	inc::Int
	stab::Float64
	allsdtol::Float64
	iStart::Float64
	iTot::Float64
	iMin::Float64
	iMax::Float64
	nImax::Int
	loadCase::LoadCase
	output::Array{Output,1}
	"""

		Step(name::AbstractString, lc::LoadCase, out::Array{Output,1})

	"""
	function Step(name::AbstractString,lc::LoadCase,out::Array{Output,1})
		new(name,false,1000,0.001,0.05,0.001,1.0,1.0e-35,0.1,20,lc,out)
	end
	"""

		Step(name::AbstractString, lc::LoadCase)

	"""
	function Step(name::AbstractString,lc::LoadCase)
		out = [Output("Node","U"),
			   Output("Element","E"),
			   Output("Element","S")]
		new(name,false,1000,0.001,0.05,0.001,1.0,1.0e-35,0.1,20,lc,out)
	end
end

"""

	show(io::IO, n::Node)

"""
function show(io::IO,n::Node)
	print(io,"Node(<$(n.num)>, <$(n.coords[1]), $(n.coords[2]), $(n.coords[3])>)")
end

"""

	show(io::IO, n::GlobNode)

"""
function show(io::IO,n::GlobNode)
	print(io,"GlobNode(<$(n.glob)@glob>, <$(n.node.num)@$(n.instance)>, <$(n.node.coords[1]), $(n.node.coords[2]), $(n.node.coords[3])>)")
end

"""

	show(io::IO, a::Array{GlobNode,1})

"""
function show(io::IO,a::Array{GlobNode,1})
	print(io,"$(length(a)) Global nodes")
end

"""

	show(io::IO, a::Array{Node,1})

"""
function show(io::IO,a::Array{Node,1})
	print(io,"$(length(a)) Nodes")
end

"""

	show(io::IO, a::Array{Instance,1})

"""
show(io::IO,a::Array{Instance,1}) = print(io,"$(length(a)) Instances")

"""

	show(io::IO, a::Array{Equation,1})

"""
show(io::IO,a::Array{Equation,1}) = print(io,"$(length(a)) Equations")

"""

	show(io::IO, a::Array{BoundCon,1})

"""
show(io::IO,a::Array{BoundCon,1}) = print(io,"$(length(a)) BoundCons")

"""

	show(io::IO, s::Step)

"""
show(io::IO,s::Step) = print(io,"Step($(s.name), $(s.loadCase.bound))")

"""

	show(io::IO, a::Array{Step,1})

"""
show(io::IO,a::Array{Step,1}) = print(io,"$(length(a)) Steps")

const float_re = "(-?\\d+.\\d*(e(\\+|-)\\d+)?)"

"""

	node_re

Defines RegEx for finding Nodes in Input file.
"""
const node_re = Regex("\\s*(\\d+),\\s*"*float_re*",\\s*"*float_re*",\\s*"*float_re)

"""

	matchIf(r::Regex, s::AbstractString, i::AbstractString)

"""
function matchIf(r::Regex,s::AbstractString,i::AbstractString)
	m = match(r,s)
	if m == nothing
		return i
	else
		return m.captures[1]
	end
end

"""

	round_node(node::Node, d::Int64)

Round the coordinates of Node node to the number of digits d.
"""
function round_node(node::Node, d::Int64)
	new_coords = round.(node.coords,digits=d)
	return Node(node.num, new_coords)
end

"""
	
	rotate(c::Vector, rot_axis::Vector, angle::Number)

"""
function rotate(c::Vector, rot_axis::Vector, angle::Number)
	angle_rad = deg2rad(angle)
	R_1 = rot_axis * rot_axis' .* (1 - cos(angle_rad))
	R_2 = [	+cos(angle_rad)				-rot_axis[3]*sin(angle_rad) +rot_axis[2]*sin(angle_rad);
		   +rot_axis[3]*sin(angle_rad)	+cos(angle_rad)				-rot_axis[1]*sin(angle_rad);
		   -rot_axis[2]*sin(angle_rad)	+rot_axis[1]*sin(angle_rad)	+cos(angle_rad)]
	R = R_1 + R_2
	return R*c
end

"""
	
rotate(c::Vector, start_rot_axis::Vector, end_rot_axis::Vector, angle::Number)

"""
function rotate(c::Vector, start_rot_axis::Vector, end_rot_axis::Vector, angle::Number)
	rot_axis = end_rot_axis - start_rot_axis
	angle_rad = deg2rad(angle)
	R_1 = rot_axis * rot_axis' .* (1 - cos(angle_rad))
	R_2 = [	+cos(angle_rad)				-rot_axis[3]*sin(angle_rad) +rot_axis[2]*sin(angle_rad);
		   +rot_axis[3]*sin(angle_rad)	+cos(angle_rad)				-rot_axis[1]*sin(angle_rad);
		   -rot_axis[2]*sin(angle_rad)	+rot_axis[1]*sin(angle_rad)	+cos(angle_rad)]
	R = R_1 + R_2
	return R * (c - start_rot_axis) + start_rot_axis
end

"""

	rotate_re

Regular Expression for loading rotational information from instances
"""

const rotate_re = Regex("\\s*"*float_re*(",\\s*"*float_re)^6)

"""

	load_rotate(rot_line::AbstractString)

"""
function load_rotate(rot_line::AbstractString)
	match_bool = [true,false,false,true,false,false,true,false,false,true,false,false,true,false,false,true,false,false,true,false,false]
	m = match(rotate_re,rot_line).captures[match_bool]
	@show m
	vals = parse.(Float64,m)
	start_rot_axis = vals[1:3]
	end_rot_axis = vals[4:6]
	angle_deg = vals[7]
	return start_rot_axis, end_rot_axis, angle_deg
end

function translate(node::Node,trans_vec::Vector)
	return Node(node.num, node.coords+trans_vec)
end

function rotate(node::Node, rot_vec::Vector, angle_deg::Number)
	new_coords = rotate(node.coords, rot_vec, angle_deg)
	return Node(node.num,new_coords)
end

function rotate(node::Node, rot_start::Vector, rot_end::Vector, angle_deg::Number)
	rot_vec = rot_end-rot_start
	return translate(rotate(translate(node,-rot_start),rot_vec,angle_deg),rot_start)
end

"""

	ToleranceError

"""
struct ToleranceError <: Exception end

"""

	NodeError

"""
struct NodeError <: Exception end

"""

	AxisError

"""
struct AxisError <: Exception end

"""

	DimensionError

"""
struct DimensionError <: Exception end

"""

	loadNode(nodeLine::Line)

"""
function loadNode(nodeLine::Line)
	props = match(node_re,value(nodeLine)).captures[[true,true,false,false,true,false,false,true,false,false]]
	return Node(parse(Int,props[1]), round.([parse(Float64,props[2]), parse(Float64,props[3]), parse(Float64,props[4])],digits=5))
end

"""

	Part(name::String, nodes::Array{Node,1})

"""
mutable struct Part
	name::String
	nodes::Array{Node,1}
end

"""

	show(io::IO, a::Dict{String,Part})

"""
function show(io::IO,a::Dict{String,Part})
	print(io,"$(length(a)) Parts")
end

"""

	findParts(f::Line)

"""
function findParts(f::File)
	return Lines(f,r"\*Part")
end

"""

	read_part(part_line::Line)

"""
function read_part(part_line::Line)
	part_name = match(r"name=(.*),?",value(part_line)).captures[1]
	nodes = Node[]
	part = File(value(Range(part_line.file, from=part_line, to=r"\*End Part")))
	node_lines = Lines(part,node_re, after=r"Node",before=r"\*")
	for nl in node_lines
		push!(nodes,loadNode(nl))
	end
	return Part(part_name, nodes)
end

"""

	read_parts(f::File)

"""
function read_parts(f::File)
	parts = Dict{String,Part}()
	p = map(findParts(f)) do l
		read_part(l)
	end
	for part in p
		parts[part.name] = part
	end
	return parts
end
"""

	findInstances(f::File)

"""
function findInstances(f::File)
	return Lines(f,r"\*Instance")
end

const trans_re = r"\s*(-?\d+.\d*(e(\+|-)\d+)?),\s*(-?\d+.\d*(e(\+|-)\d+)?),\s*(-?\d+.\d*(e(\+|-)\d+)?)"

"""

	readInstance(instLine::Line)

"""
function readInstance(instLine::Line)
	inst = Range(instLine.file,from=instLine,to=r"\*End Instance")
	name = match(r"Instance, name=(.+), part=.+",value(instLine)).captures[1]
	@info "Reading Instance $(name)"
	part = match(r"Instance, name=.+, part=(.+)",value(instLine)).captures[1]
	if length(inst) == 2
		return Instance(part, name)
	elseif length(inst) == 3
		trans = parse.(Float64,match(trans_re,value(inst[2])).captures[[true,false,false,true,false,false,true,false,false]])
		return Instance(part, name, trans)
	elseif length(inst) == 4
		trans = parse.(Float64,match(trans_re,value(inst[2])).captures[[true,false,false,true,false,false,true,false,false]])
		rot_start, rot_end, angle_deg = load_rotate(value(inst[3]))
		return Instance(part, name, trans, rot_start, rot_end, angle_deg)
	else
		@error "Instance $(name) has $(length(inst)) lines!"
	end
end

"""

	readInstances(f::File)

"""
function readInstances(f::File)
	i = map(findInstances(f)) do l
		readInstance(l)
	end
	return i
end

"""

	loadGlobNodes(f::File)

"""
function loadGlobNodes(f::File)
	parts = read_parts(f)
	instances = readInstances(f)
	globCounter = 1
	globNodes = Array{GlobNode,1}()
	for i in instances
		part_nodes = parts[i.part].nodes
		inst_nodes = map(part_nodes) do p_n
			round_node(rotate(translate(p_n,i.trans),i.rot_start,i.rot_end,i.angle_deg),4)
		end
		for i_n in inst_nodes
			push!(globNodes, GlobNode(globCounter,i.name,i_n))
			globCounter += 1
		end
		i.nodes = inst_nodes
	end
	return parts, instances, globNodes
end

"""

	ContactPair(master::String, slave::String)

"""
mutable struct ContactPair
	name::String
	master::String
	slave::String
end

show(io::IO, cp::ContactPair) = print("ContactPair($(cp.name) Master:$(cp.master) <-> Slave:$(cp.slave))")

"""

	find_ties(f::File)

"""
function find_ties(f::File)
	return Lines(f,r"\*Tie")
end

"""

	read_tie(part_line::Line)

"""
function read_tie(tie_line::Line)
	tie_name = match(r"name=(.*),",value(tie_line)).captures[1]
	surfaces = match(r"(.+),\s+(.+)",value(tie_line+1)).captures
	return ContactPair(tie_name, surfaces[1], surfaces[2])
end

"""

	read_ties(f::File)

"""
function read_ties(f::File)
	ties = Dict{String,ContactPair}()
	t = map(find_ties(f)) do l
		read_tie(l)
	end
	for tie in t
		ties[tie.name] = tie
	end
	return ties
end

"""

	Surface(name::String, sets::Array{Tuple{String, String}})

"""
mutable struct Surface
	name::String
	sets::Array{Tuple{String, String},1}
end

const surf_re = r"(.*),\s(.*)"

"""

	find_surfaces(f::File)

"""
function find_surfaces(f::File)
	return Lines(f,r"\*Surface,")
end

"""

	read_surface(surf_line::Line)

"""
function read_surface(surf_line::Line)
	name = match(r"name=(.*)",value(surf_line)).captures[1]
	sets = Array{Tuple{String,String},1}()
	surf_def = Lines(surf_line.file,surf_re,after=surf_line, before=Line(surf_line.file,r"\*",after=surf_line))
	sets = map(surf_def) do l
		def = match(surf_re,value(l)).captures
		(def[1],def[2])
	end
	return Surface(name, sets)
end

"""

	read_surfaces(f::File)

"""
function read_surfaces(f::File)
	surfaces = Dict{String,Surface}()
	s = map(find_surfaces(f)) do l
		read_surface(l)
	end
	for surface in s
		surfaces[surface.name] = surface
	end
	return surfaces
end

const node_num = Dict("C3D8" => 8)

const elem_re = Dict("C3D8" => r"\s*(\d+),\s+(\d+),\s+(\d+),\s+(\d+),\s+(\d+),\s+(\d+),\s+(\d+),\s+(\d+),\s+(\d+)")

"""

	ElementType(type::String, nn::Int64)

"""
struct ElementType
	type::String
	nn::Int64
	"""

		ElementType(type::name)

	"""
	function ElementType(type::String)
		return new(type, node_num[type])
	end
end


"""

	Element(type::String, num::Int64, nodes::Array{Node,1})

"""
mutable struct Element
	type::ElementType
	num::Int64
	nodes::Array{Node,1}
end

"""

	read_element(elem_line::Line, type::String)

"""
function read_element(elem_line::Line, type::String)
	props = match(elem_re[type],value(elem_line)).captures
	nodes = props[2:end]
	return Element(type, parse(Int,props[1]), parse.(Int64,nodes))
end


"""

	ElementSet(name::String, elements::Array{Element,1})

"""
mutable struct ElementSet
	name::String
	elements::Array{Element,1}
end

"""

	find_elsets(f::File)

"""
function find_elsets(f::File)
	return Lines(f,r"\*Surface,")
end

"""

	read_elset(surf_line::Line)

"""
function read_elset(surf_line::Line)
	name = match(r"name=(.*)",value(surf_line)).captures[1]
	sets = Array{Tuple{String,String},1}()
	surf_def = Lines(surf_line.file,surf_re,after=surf_line, before=Line(surf_line.file,r"\*",after=surf_line))
	sets = map(surf_def) do l
		def = match(surf_re,value(l)).captures
		(def[1],def[2])
	end
	return Surface(name, sets)
end

"""

	read_elsets(f::File)

"""
function read_elsets(f::File)
	elsets = Dict{String,Surface}()
	e = map(find_elsets(f)) do l
		read_elsets(l)
	end
	for elset in e
		elsets[elset.name] = elset
	end
	return elsets
end

function count_char(str::String,char::Char)
	counter = 0
	for c in str
		if c == char
			counter +=1
		end
	end
	return counter
end

const slave_re_1 = "\\s*(\\d+)"
const slave_re_2 = ",\\s*(\\d+)"

function collect_slaves(f)
	try
		set_lines = Lines(f,r"\*Nset, nset=Slaves, instance")
		slaves = Dict{String,Array{Int64,1}}()
		for sl in set_lines
			instance = match(r"instance=(.*),?",value(sl)).captures[1]
			node_def = Lines(f, r"\s*\d+,", after=sl, before=Line(f,r"\*",after=sl))
			slave_nodes = Array{Int64,1}()
			for nl in node_def
				n = count_char(value(nl),',')
				re = Regex(slave_re_1*slave_re_2^n)
				nodes = parse.(Int64,match(re,value(nl)).captures)
				append!(slave_nodes,nodes)
			end
			slaves[instance] = slave_nodes
		end
		return slaves
	catch
		return Dict("Slaves"=>Int64[])
	end
end

"""

	getLength(nodes::Array{GlobNode})

"""
function getLength(nodes::Array{GlobNode})
	x = map(nodes) do n
		n.node.coords[1]
	end
	y = map(nodes) do n
		n.node.coords[2]
	end
	z = map(nodes) do n
		n.node.coords[3]
	end
	min = [minimum(x), minimum(y), minimum(z)]
	max = [maximum(x), maximum(y), maximum(z)]
	return min, max, max-min
end

include("abqModel.jl")
include("nodeDesignation.jl")
include("periodicBoundaryConditions.jl")
include("pbc2D.jl")
include("pbc3D.jl")
include("loadBoundaryConditions.jl")

end # module
