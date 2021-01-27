"""

	pbc1d!(abq::AbqModel)

Generate periodic boundary conditions for onedimensional periodicity and
append them to the given AbqModel.
"""
function pbc1D!(abq::AbqModel)
	# abbreviation for coordinate system
	c = abq.csys
	# abbreviation for unit cell dimensions
	l = abq.dim
	# initiate counter for equation numbering
	i=1
	################
	# VERTEX NODES #
	################
	# no equations for vertices needed!
	##############
	# EDGE NODES #
	##############
	for n = 1:length(abq.edges["ST"])
		this_node = abq.edges["ST"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		x = this_node.node.coords
		for a = 1:3
			i+=1
			push!(abq.eqns,Equation(i,["ST-$(n)","SB-$(n)","SWT","SWB","SET","SEB"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,x[c[2]]/l[c[2]]-1.0,1.0-x[c[2]]/l[c[2]],-x[c[2]]/l[c[2]],x[c[2]]/l[c[2]]]))
		end
	end
	for n = 1:length(abq.edges["WT"])
		this_node = abq.edges["WT"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		x = this_node.node.coords
		for a = 1:3
			i+=1
			push!(abq.eqns,Equation(i,["WT-$(n)","WB-$(n)","SWT","SWB","NWT","NWB"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,x[c[3]]/l[c[3]]-1.0,1.0-x[c[3]]/l[c[3]],-x[c[3]]/l[c[3]],x[c[3]]/l[c[3]]]))
		end
	end
	for n = 1:length(abq.edges["ET"])
		this_node = abq.edges["ET"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		x = this_node.node.coords
		for a = 1:3
			i+=1
			push!(abq.eqns,Equation(i,["ET-$(n)","EB-$(n)","SET","SEB","NET","NEB"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,x[c[3]]/l[c[3]]-1.0,1.0-x[c[3]]/l[c[3]],-x[c[3]]/l[c[3]],x[c[3]]/l[c[3]]]))
		end
	end
	for n = 1:length(abq.edges["NT"])
		this_node = abq.edges["NT"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		x = this_node.node.coords
		for a = 1:3
			i+=1
			push!(abq.eqns,Equation(i,["NT-$(n)","NB-$(n)","NWT","NWB","NET","NEB"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,x[c[2]]/l[c[2]]-1.0,1.0-x[c[2]]/l[c[2]],-x[c[2]]/l[c[2]],x[c[2]]/l[c[2]]]))
		end
	end
	##############
	# FACE NODES #
	##############
	for n = 1:length(abq.faces["T"])
		this_node = abq.faces["T"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		x = this_node.node.coords
		for a = 1:3
			i+=1
			push!(abq.eqns,Equation(i,["T-$(n)","B-$(n)","SWT","SWB","SET","SEB","NWT","NWB","NET","NEB"],
									  [c[a],c[a],c[a],c[a],c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,-1.0+x[c[2]]/l[c[2]]+x[c[2]]/l[c[2]]-x[c[2]]/l[c[2]]*x[c[3]]/l[c[3]],1.0-x[c[2]]/l[c[2]]-x[c[2]]/l[c[2]]+x[c[2]]/l[c[2]]*x[c[3]]/l[c[3]],
									   -x[c[2]]/l[c[2]]+x[c[2]]/l[c[2]]*x[c[3]]/l[c[3]],x[c[2]]/l[c[2]]-x[c[2]]/l[c[2]]*x[c[3]]/l[c[3]],
									   -x[c[3]]/l[c[3]]+x[c[2]]/l[c[2]]*x[c[3]]/l[c[3]],x[c[3]]/l[c[3]]-x[c[2]]/l[c[2]]*x[c[3]]/l[c[3]],
									   -x[c[2]]/l[c[2]]*x[c[3]]/l[c[3]],x[c[2]]/l[c[2]]*x[c[3]]/l[c[3]]]))
		end
	end
	@info "$(i) equations written to AbqModel."
	return
end
