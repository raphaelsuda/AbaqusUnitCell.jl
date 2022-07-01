"""

	pbc2d!(abq::AbqModel)

Generate periodic boundary conditions for twodimensional periodicity and
append them to the given AbqModel.
"""
function pbc2D!(abq::AbqModel)
	# abbreviation for coordinate system
	c = abq.csys
	# abbreviation for unit cell dimensions
	l = abq.dim
	# initiate counter for equation numbering
	i=1
	##
	## VERTICES
	##
	push!(abq.eqns,Equation(i,["NWT","SWT","NWB","SWB"],[c[1],c[1],c[1],c[1]],[1.0,-1.0,-1.0,1.0]))
	i+=1
	push!(abq.eqns,Equation(i,["SET","SWT","SEB","SWB"],[c[1],c[1],c[1],c[1]],[1.0,-1.0,-1.0,1.0]))
	i+=1
	push!(abq.eqns,Equation(i,["NEB","NWB","SEB","SWB","SWB","SWT","SEB","SET","SWB","SWT","NWB","NWT"],
							[c[1],c[1],c[1],c[1],c[2],c[2],c[2],c[2],c[3],c[3],c[3],c[3]],
							[1.0,-1.0,-1.0,1.0,
							 l[c[2]]/l[c[1]],-l[c[2]]/l[c[1]],-l[c[2]]/l[c[1]],l[c[2]]/l[c[1]],
							 l[c[3]]/l[c[1]],-l[c[3]]/l[c[1]],-l[c[3]]/l[c[1]],l[c[3]]/l[c[1]]]))
	for a = 2:3
		i+=1
		push!(abq.eqns,Equation(i,["NEB","NWB","SEB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
	end	
	i+=1
	push!(abq.eqns,Equation(i,["NET","NWT","SEB","SWB","SWB","SWT","SEB","SET","SWB","SWT","NWB","NWT"],
							[c[1],c[1],c[1],c[1],c[2],c[2],c[2],c[2],c[3],c[3],c[3],c[3]],
							[1.0,-1.0,-1.0,1.0,
							 l[c[2]]/l[c[1]],-l[c[2]]/l[c[1]],-l[c[2]]/l[c[1]],l[c[2]]/l[c[1]],
							 l[c[3]]/l[c[1]],-l[c[3]]/l[c[1]],-l[c[3]]/l[c[1]],l[c[3]]/l[c[1]]]))
	for a = 2:3
		i+=1
		push!(abq.eqns,Equation(i,["NET","NWT","SET","SWT"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
	end	
	##
	## EDGES
	##
	for n = 1:length(abq.edges["SE"])
		this_node = abq.edges["SE"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		i+=1
		x = abq.edges["SE"][n].node.coords
		push!(abq.eqns,Equation(i,["SE-$(n)","SW-$(n)","SEB","SWB"],[c[1],c[1],c[1],c[1]],[1.0,-1.0,-1.0,1.0]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["SE-$(n)","SW-$(n)","SEB","SWB","SWT","SET"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,(x[c[1]]-abq.minC[c[1]])/l[c[1]]-1.0,1.0-(x[c[1]]-abq.minC[c[1]])/l[c[1]],(x[c[1]]-abq.minC[c[1]])/l[c[1]],-(x[c[1]]-abq.minC[c[1]])/l[c[1]]]))
		end
	end
	for n = 1:length(abq.edges["NW"])
		this_node = abq.edges["NW"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		i+=1
		x = abq.edges["NW"][n].node.coords
		push!(abq.eqns,Equation(i,["NW-$(n)","SW-$(n)","NWB","SWB"],[c[1],c[1],c[1],c[1]],[1.0,-1.0,-1.0,1.0]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["NW-$(n)","SW-$(n)","NWB","SWB","SWT","NWT"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,(x[c[1]]-abq.minC[c[1]])/l[c[1]]-1.0,1.0-(x[c[1]]-abq.minC[c[1]])/l[c[1]],(x[c[1]]-abq.minC[c[1]])/l[c[1]],-(x[c[1]]-abq.minC[c[1]])/l[c[1]]]))
		end
	end
	for n = 1:length(abq.edges["NE"])
		this_node = abq.edges["NE"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		i+=1
		x = abq.edges["NE"][n].node.coords
		push!(abq.eqns,Equation(i,["NE-$(n)","NW-$(n)","SEB","SWB","SWB","SWT","SEB","SET","SWB","SWT","NWB","NWT"],
								  [c[1],c[1],c[1],c[1],c[2],c[2],c[2],c[2],c[3],c[3],c[3],c[3]],
								  [1.0,-1.0,-1.0,1.0,
								   l[c[2]]/l[c[1]],-l[c[2]]/l[c[1]],-l[c[2]]/l[c[1]],l[c[2]]/l[c[1]],
								   l[c[3]]/l[c[1]],-l[c[3]]/l[c[1]],-l[c[3]]/l[c[1]],l[c[3]]/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["NE-$(n)","NW-$(n)","SEB","SWB","SWT","SET"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,(x[c[1]]-abq.minC[c[1]])/l[c[1]]-1.0,1.0-(x[c[1]]-abq.minC[c[1]])/l[c[1]],(x[c[1]]-abq.minC[c[1]])/l[c[1]],-(x[c[1]]-abq.minC[c[1]])/l[c[1]]]))
		end
	end
	for n = 1:length(abq.edges["EB"])
		this_node = abq.edges["EB"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		i+=1
		x = abq.edges["EB"][n].node.coords
		push!(abq.eqns,Equation(i,["EB-$(n)","WB-$(n)","SEB","SWB","SWB","SWT","SEB","SET"],
								  [c[1],c[1],c[1],c[1],c[2],c[2],c[2],c[2]],
								  [1.0,-1.0,-1.0,1.0,(x[c[2]]-abq.minC[c[2]])/l[c[1]],-(x[c[2]]-abq.minC[c[2]])/l[c[1]],-(x[c[2]]-abq.minC[c[2]])/l[c[1]],(x[c[2]]-abq.minC[c[2]])/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["EB-$(n)","WB-$(n)","SEB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
		end
	end
	for n = 1:length(abq.edges["ET"])
		this_node = abq.edges["ET"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		i+=1
		x = abq.edges["ET"][n].node.coords
		push!(abq.eqns,Equation(i,["ET-$(n)","WT-$(n)","SEB","SWB","SWB","SWT","SEB","SET"],
								  [c[1],c[1],c[1],c[1],c[2],c[2],c[2],c[2]],
								  [1.0,-1.0,-1.0,1.0,(x[c[2]]-abq.minC[c[2]])/l[c[1]],-(x[c[2]]-abq.minC[c[2]])/l[c[1]],-(x[c[2]]-abq.minC[c[2]])/l[c[1]],(x[c[2]]-abq.minC[c[2]])/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["ET-$(n)","WT-$(n)","SET","SWT"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
		end
	end
	for n = 1:length(abq.edges["NB"])
		this_node = abq.edges["NB"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		i+=1
		x = abq.edges["NB"][n].node.coords
		push!(abq.eqns,Equation(i,["NB-$(n)","SB-$(n)","NWB","SWB","SWB","SWT","NWB","NWT"],
								  [c[1],c[1],c[1],c[1],c[3],c[3],c[3],c[3]],
								  [1.0,-1.0,-1.0,1.0,(x[c[3]]-abq.minC[c[3]])/l[c[1]],-(x[c[3]]-abq.minC[c[3]])/l[c[1]],-(x[c[3]]-abq.minC[c[3]])/l[c[1]],(x[c[3]]-abq.minC[c[3]])/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["NB-$(n)","SB-$(n)","NWB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
		end
	end
	for n = 1:length(abq.edges["NT"])
		this_node = abq.edges["NT"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		i+=1
		x = abq.edges["NT"][n].node.coords
		push!(abq.eqns,Equation(i,["NT-$(n)","ST-$(n)","NWB","SWB","SWB","SWT","NWB","NWT"],
								  [c[1],c[1],c[1],c[1],c[3],c[3],c[3],c[3]],
								  [1.0,-1.0,-1.0,1.0,(x[c[3]]-abq.minC[c[3]])/l[c[1]],-(x[c[3]]-abq.minC[c[3]])/l[c[1]],-(x[c[3]]-abq.minC[c[3]])/l[c[1]],(x[c[3]]-abq.minC[c[3]])/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["NT-$(n)","ST-$(n)","NWT","SWT"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
		end
	end
	##
	## SURFACES
	##
	for n = 1:length(abq.faces["E"])
		this_node = abq.faces["E"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		i+=1
		x = abq.faces["E"][n].node.coords
		push!(abq.eqns,Equation(i,["E-$(n)","W-$(n)","SEB","SWB","SWB","SWT","SEB","SET"],
								  [c[1],c[1],c[1],c[1],c[2],c[2],c[2],c[2]],
								  [1.0,-1.0,-1.0,1.0,(x[c[2]]-abq.minC[c[2]])/l[c[1]],-(x[c[2]]-abq.minC[c[2]])/l[c[1]],-(x[c[2]]-abq.minC[c[2]])/l[c[1]],(x[c[2]]-abq.minC[c[2]])/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["E-$(n)","W-$(n)","SEB","SWB","SWT","SET"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,(x[c[1]]-abq.minC[c[1]])/l[c[1]]-1.0,1.0-(x[c[1]]-abq.minC[c[1]])/l[c[1]],(x[c[1]]-abq.minC[c[1]])/l[c[1]],-(x[c[1]]-abq.minC[c[1]])/l[c[1]]]))
		end
	end
	for n = 1:length(abq.faces["N"])
		this_node = abq.faces["N"][n]
		if this_node.instance in keys(abq.slaves)
			if this_node.node.num in abq.slaves[this_node.instance]
				continue
			end
		end
		i+=1
		x = abq.faces["N"][n].node.coords
		push!(abq.eqns,Equation(i,["N-$(n)","S-$(n)","NWB","SWB","SWB","SWT","NWB","NWT"],
								  [c[1],c[1],c[1],c[1],c[3],c[3],c[3],c[3]],
								  [1.0,-1.0,-1.0,1.0,(x[c[3]]-abq.minC[c[3]])/l[c[1]],-(x[c[3]]-abq.minC[c[3]])/l[c[1]],-(x[c[3]]-abq.minC[c[3]])/l[c[1]],(x[c[3]]-abq.minC[c[3]])/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["N-$(n)","S-$(n)","NWB","SWB","SWT","NWT"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,(x[c[1]]-abq.minC[c[1]])/l[c[1]]-1.0,1.0-(x[c[1]]-abq.minC[c[1]])/l[c[1]],(x[c[1]]-abq.minC[c[1]])/l[c[1]],-(x[c[1]]-abq.minC[c[1]])/l[c[1]]]))
		end
	end
	@info "$(i) equations written to AbqModel."
	return
end
