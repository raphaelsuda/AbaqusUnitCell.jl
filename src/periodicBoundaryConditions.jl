"""

	generate(eq::Equation)

"""
function generate(eq::Equation)
	n = length(eq.nodes)
	eqString = ["** Constraint: Eqn-$(eq.num)",
				"*Equation",
				"$(n)"]
	for i = 1:n
		push!(eqString,"$(eq.nodes[i]), $(eq.dof[i]), $(eq.coef[i])")
	end
	return eqString
end

"""

	pbc!(abq::AbqModel)

"""
function pbc!(abq::AbqModel)
	abq.defDim ? println("Periodicity not explicitly defined. Default value $(abq.pbcdim)-dimensional periodicity is used.") : print("")
	if abq.pbcdim == 1
		pbc1D!(abq)
	elseif abq.pbcdim == 2
		pbc2D!(abq)
	elseif abq.pbcdim == 3
		pbc3D!(abq)
	end
return
end
