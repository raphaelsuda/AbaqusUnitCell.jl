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
	# TODO equations for edge nodes
	for n = 1:length(abq.edges["SE"])
		i+=1
		x = abq.edges["SE"][n].coords
		push!(abq.eqns,Equation(i,["SE-$(n)","SW-$(n)","SEB","SWB"],[c[1],c[1],c[1],c[1]],[1.0,-1.0,-1.0,1.0]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["SE-$(n)","SW-$(n)","SEB","SWB","SWT","SET"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,x[c[1]]/l[c[1]]-1.0,1.0-x[c[1]]/l[c[1]],x[c[1]]/l[c[1]],-x[c[1]]/l[c[1]]]))
		end
	end
	for n = 1:length(abq.edges["NW"])
		i+=1
		x = abq.edges["NW"][n].coords
		push!(abq.eqns,Equation(i,["NW-$(n)","SW-$(n)","NWB","SWB"],[c[1],c[1],c[1],c[1]],[1.0,-1.0,-1.0,1.0]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["NW-$(n)","SW-$(n)","NWB","SWB","SWT","NWT"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,x[c[1]]/l[c[1]]-1.0,1.0-x[c[1]]/l[c[1]],x[c[1]]/l[c[1]],-x[c[1]]/l[c[1]]]))
		end
	end
	for n = 1:length(abq.edges["NE"])
		i+=1
		x = abq.edges["NE"][n].coords
		push!(abq.eqns,Equation(i,["NE-$(n)","NW-$(n)","SEB","SWB","SWB","SWT","SEB","SET","SWB","SWT","NWB","NWT"],
								  [c[1],c[1],c[1],c[1],c[2],c[2],c[2],c[2],c[3],c[3],c[3],c[3]],
								  [1.0,-1.0,-1.0,1.0,
								   l[c[2]]/l[c[1]],l[c[2]]/l[c[1]],l[c[2]]/l[c[1]],l[c[2]]/l[c[1]],
								   l[c[3]]/l[c[1]],l[c[3]]/l[c[1]],l[c[3]]/l[c[1]],l[c[3]]/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["NE-$(n)","NW-$(n)","SEB","SWB","SWT","SET"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,x[c[1]]/l[c[1]]-1.0,1.0-x[c[1]]/l[c[1]],x[c[1]]/l[c[1]],-x[c[1]]/l[c[1]]]))
		end
	end
	for n = 1:length(abq.edges["EB"])
		i+=1
		x = abq.edges["EB"][n].coords
		push!(abq.eqns,Equation(i,["EB-$(n)","WB-$(n)","SEB","SWB","SWB","SWT","SEB","SET"],
								  [c[1],c[1],c[1],c[1],c[2],c[2],c[2],c[2]],
								  [1.0,-1.0,-1.0,1.0,x[c[2]]/l[c[1]],-x[c[2]]/l[c[1]],-x[c[2]]/l[c[1]],x[c[2]]/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["EB-$(n)","WB-$(n)","SEB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
		end
	end
	for n = 1:length(abq.edges["ET"])
		i+=1
		x = abq.edges["ET"][n].coords
		push!(abq.eqns,Equation(i,["ET-$(n)","WT-$(n)","SEB","SWB","SWB","SWT","SEB","SET"],
								  [c[1],c[1],c[1],c[1],c[2],c[2],c[2],c[2]],
								  [1.0,-1.0,-1.0,1.0,x[c[2]]/l[c[1]],-x[c[2]]/l[c[1]],-x[c[2]]/l[c[1]],x[c[2]]/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["ET-$(n)","WT-$(n)","SEB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
		end
	end
	for n = 1:length(abq.edges["NB"])
		i+=1
		x = abq.edges["NB"][n].coords
		push!(abq.eqns,Equation(i,["NB-$(n)","SB-$(n)","NWB","SWB","SWB","SWT","NWB","NWT"],
								  [c[1],c[1],c[1],c[1],c[3],c[3],c[3],c[3]],
								  [1.0,-1.0,-1.0,1.0,x[c[3]]/l[c[1]],-x[c[3]]/l[c[1]],-x[c[3]]/l[c[1]],x[c[3]]/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["NB-$(n)","SB-$(n)","NWB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
		end
	end
	for n = 1:length(abq.edges["NT"])
		i+=1
		x = abq.edges["NT"][n].coords
		push!(abq.eqns,Equation(i,["NT-$(n)","ST-$(n)","NWB","SWB","SWB","SWT","NWB","NWT"],
								  [c[1],c[1],c[1],c[1],c[3],c[3],c[3],c[3]],
								  [1.0,-1.0,-1.0,1.0,x[c[3]]/l[c[1]],-x[c[3]]/l[c[1]],-x[c[3]]/l[c[1]],x[c[3]]/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["NT-$(n)","ST-$(n)","NWB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
		end
	end
	for n = 1:length(abq.faces["E"])
		i+=1
		x = abq.faces["E"][n].coords
		push!(abq.eqns,Equation(i,["E-$(n)","W-$(n)","SEB","SWB","SWB","SWT","SEB","SET"],
								  [c[1],c[1],c[1],c[1],c[2],c[2],c[2],c[2]],
								  [1.0,-1.0,-1.0,1.0,x[c[2]]/l[c[1]],-x[c[2]]/l[c[1]],-x[c[2]]/l[c[1]],x[c[2]]/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["E-$(n)","W-$(n)","SEB","SWB","SWT","SET"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,x[c[1]]/l[c[1]]-1.0,1.0-x[c[1]]/l[c[1]],x[c[1]]/l[c[1]],-x[c[1]]/l[c[1]]]))
		end
	end
	for n = 1:length(abq.faces["N"])
		i+=1
		x = abq.faces["N"][n].coords
		push!(abq.eqns,Equation(i,["N-$(n)","S-$(n)","NWB","SWB","SWB","SWT","NWB","NWT"],
								  [c[1],c[1],c[1],c[1],c[3],c[3],c[3],c[3]],
								  [1.0,-1.0,-1.0,1.0,x[c[3]]/l[c[1]],-x[c[3]]/l[c[1]],-x[c[3]]/l[c[1]],x[c[3]]/l[c[1]]]))
		for a = 2:3
			i+=1
			push!(abq.eqns,Equation(i,["N-$(n)","S-$(n)","NWB","SWB","SWT","NWT"],
									  [c[a],c[a],c[a],c[a],c[a],c[a]],
									  [1.0,-1.0,x[c[1]]/l[c[1]]-1.0,1.0-x[c[1]]/l[c[1]],x[c[1]]/l[c[1]],-x[c[1]]/l[c[1]]]))
		end
	end
	println("$(i) equations written to AbqModel.")
	return
end
