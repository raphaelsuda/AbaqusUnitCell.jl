"""
	pbc3D!(abq::AbqModel)

	Generate periodic boundary conditions for threedimensional periodicity and
	append them to the given AbqModel.
"""
function pbc3D!(abq::AbqModel)
	# abbreviation for coordinate system
	c = abq.csys
	# abbreviation for unit cell dimensions
	l = abq.dim
	# initiate counter for equation numbering
	i=1
	######################
	## vertex equations ##
	######################
	
	# vertex NEB
	for a = 1:3
		push!(abq.eqns,Equation(i,["NEB","NWB","SEB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
		i+=1
	end
	# vertex SET
	for a = 1:3
		push!(abq.eqns,Equation(i,["SET","SWT","SEB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
		i+=1
	end
	# vertex NWT
	for a = 1:3
		push!(abq.eqns,Equation(i,["NWT","SWT","NWB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
		i+=1
	end
	# vertex NET
	for a = 1:3
		push!(abq.eqns,Equation(i,["NET","SWT","NWB","SEB","SWB"],[c[a],c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,-1.0,2.0]))
		i+=1
	end
	
	####################
	## edge equations ##
	####################
	
	# edge NB
	for n = 1:length(abq.edges["NB"])
		x = abq.edges["NB"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["NB-$(n)","SB-$(n)","NWB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
			i+=1
		end
	end
	# edge ST
	for n = 1:length(abq.edges["ST"])
		x = abq.edges["ST"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["ST-$(n)","SB-$(n)","SWT","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
			i+=1
		end
	end
	# edge NT
	for n = 1:length(abq.edges["NT"])
		x = abq.edges["NT"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["NT-$(n)","SB-$(n)","NWB","SWT","SWB"],[c[a],c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,-1.0,2.0]))
			i+=1
		end
	end
	# edge EB
	for n = 1:length(abq.edges["EB"])
		x = abq.edges["EB"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["EB-$(n)","WB-$(n)","SEB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
			i+=1
		end
	end
	# edge WT
	for n = 1:length(abq.edges["WT"])
		x = abq.edges["WT"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["WT-$(n)","WB-$(n)","SWT","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
			i+=1
		end
	end
	# edge ET
	for n = 1:length(abq.edges["ET"])
		x = abq.edges["ET"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["ET-$(n)","WB-$(n)","SEB","SWT","SWB"],[c[a],c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,-1.0,2.0]))
			i+=1
		end
	end
	# edge SE
	for n = 1:length(abq.edges["SE"])
		x = abq.edges["SE"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["SE-$(n)","SW-$(n)","SEB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
			i+=1
		end
	end
	# edge NW
	for n = 1:length(abq.edges["NW"])
		x = abq.edges["NW"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["NW-$(n)","SW-$(n)","NWB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
			i+=1
		end
	end
	# edge NE
	for n = 1:length(abq.edges["NE"])
		x = abq.edges["NE"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["NE-$(n)","SW-$(n)","SEB","NWB","SWB"],[c[a],c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,-1.0,2.0]))
			i+=1
		end
	end

	####################
	## face equations ##
	####################
	
	# face E
	for n = 1:length(abq.faces["E"])
		x = abq.faces["E"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["E-$(n)","W-$(n)","SEB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
			i+=1
		end
	end
	# face N
	for n = 1:length(abq.faces["N"])
		x = abq.faces["N"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["N-$(n)","S-$(n)","NWB","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
			i+=1
		end
	end
	# face T
	for n = 1:length(abq.faces["T"])
		x = abq.faces["T"][n].node.coords
		for a = 1:3
			push!(abq.eqns,Equation(i,["T-$(n)","B-$(n)","SWT","SWB"],[c[a],c[a],c[a],c[a]],[1.0,-1.0,-1.0,1.0]))
			i+=1
		end
	end
	@info "$(i) equations written to AbqModel."
	return
end
