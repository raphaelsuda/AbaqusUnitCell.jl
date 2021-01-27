"""

	lbcIniDict

"""
const lbcIniDict = Dict("SWB"=>[1,2,3],
					"SWT"=>[2,3],
					"NWB"=>[1,2,3],
					"NWT"=>[2,3],
					"SEB"=>[1,2,3],
					"SET"=>[2,3])

"""

	lbcIni(abq::AbqModel)

"""
function lbcIni(abq::AbqModel)
	bc = Array{BoundCon,1}()
	i=0
	for v in keys(lbcIniDict)
		for d in lbcIniDict[v]
			i+=1
			push!(bc,BoundCon("Disp-BC-$(i)",v,abq.csys[d]))
		end
	end
	return bc
end

"""

	generate(bc::BoundCon)

"""
function generate(bc::BoundCon)
	bcString = ["** Name: $(bc.name) Type: $(bc.type)",
				"*Boundary",
				"$(bc.node), $(bc.dof), $(bc.dof), $(bc.disp)"]
	if bc.new
		bcString[2]=string(bcString[2],", op=NEW")
	end
	return bcString
end

"""

	loadCases

"""
const loadCases = Dict("eps11"=>Dict("SWB"=>[[1,2,3],[0,0,0]],
									 "SWT"=>[[1,2,3],[1,0,0]],
									 "SEB"=>[[1,2,3],[0,0,0]],
									 "SET"=>[[1,2,3],[1,0,0]],
									 "NWB"=>[[1,2,3],[0,0,0]],
									 "NWT"=>[[1,2,3],[1,0,0]],
									 "NEB"=>[[1],[0]],
									 "NET"=>[[1],[1]]),
					   "eps22"=>Dict("SWB"=>[[1,2,3],[0,0,0]],
									 "SWT"=>[[2,3],[0,0]],
									 "SEB"=>[[1,2,3],[0,0,0]],
									 "SET"=>[[2,3],[0,0]],
									 "NWB"=>[[1,2,3],[0,2,0]],
									 "NWT"=>[[2,3],[2,0]]),
					   "eps33"=>Dict("SWB"=>[[1,2,3],[0,0,0]],
									 "SWT"=>[[2,3],[0,0]],
									 "SEB"=>[[1,2,3],[0,0,3]],
									 "SET"=>[[2,3],[0,3]],
									 "NWB"=>[[1,2,3],[0,0,0]],
									 "NWT"=>[[2,3],[0,0]]),
					   "eps12"=>Dict("SWB"=>[[1,2,3],[0,0,0]],
									 "SWT"=>[[2,3],[1,0]],
									 "SEB"=>[[1,2,3],[0,0,0]],
									 "SET"=>[[2,3],[1,0]],
									 "NWB"=>[[1,2,3],[2,0,0]],
									 "NWT"=>[[2,3],[1,0]]),
					   "eps13"=>Dict("SWB"=>[[1,2,3],[0,0,0]],
									 "SWT"=>[[2,3],[0,1]],
									 "SEB"=>[[1,2,3],[3,0,0]],
									 "SET"=>[[2,3],[0,1]],
									 "NWB"=>[[1,2,3],[0,0,0]],
									 "NWT"=>[[2,3],[0,1]]),
					   "eps23"=>Dict("SWB"=>[[1,2,3],[0,0,0]],
									 "SWT"=>[[2,3],[0,0]],
									 "SEB"=>[[1,2,3],[0,3,0]],
									 "SET"=>[[2,3],[3,0]],
									 "NWB"=>[[1,2,3],[0,0,2]],
									 "NWT"=>[[2,3],[0,2]]))

const lcAddition = Dict("eps11"=>1,"eps12"=>2,"eps13"=>3)

"""

	generate(lc::loadCase)

"""
function generate(lc::LoadCase)
	lcString = ["**","** BOUNDARY CONDITIONS","**"]
	for b in lc.bound
		append!(lcString,generate(b))
	end
	return lcString
end

"""

	generate(out::Output)

"""
function generate(out::Output)
	return ["*$(out.type) Output$(out.type=="Element" ? ", directions=YES" : "")",
			"$(out.sign)"]
end

"""

	generate(step::Step)

"""
function generate(step::Step)
# 	stepString = ["** --------------------------------------------------",
# 				  "**","** STEP: $(step.name)","**",
# 				  "*Step, name=$(step.name), nlgeom=$(step.nl ? "YES" : "NO"), inc=$(step.inc)",
# 				  "*Static, $(step.stab!=0 ? "stabilize=$(step.stab), " : "")$(step.allsdtol!=0 ? "allsdtol=$(step.allsdtol), " : "")continue=NO",
# 				  "$(step.iStart), $(step.iTot), $(step.iMin), $(step.iMax)"]
# 	append!(stepString,["**","** CONTROLS","**",
# 						"*Controls, reset",
# 						"*Controls, parameters=time incrementation",
# 						", , , , , , , $(step.nImax), , , ",
# 						", , , , , , 4.,",
# 						", 4., , , , , ,"])
# 	append!(stepString,generate(step.loadCase))
	stepString = ["** --------------------------------------------------",
				  "**","** STEP: $(step.name)","**",
				  "*Step, name=$(step.name), nlgeom=NO",
				  "*Static",
				  "$(step.iStart), $(step.iTot), $(step.iMin), $(step.iMax)"]
	append!(stepString,generate(step.loadCase))
	append!(stepString,["**","** OUTPUT REQUESTS","**"])
# 	append!(stepString,["*Output, field"])
	append!(stepString,["*Output, field, variable=PRESELECT"])
# 	for o in step.output
# 		append!(stepString,generate(o))
# 	end
 	append!(stepString, ["*Output, field"])
	for v in keys(vertices)
		append!(stepString,["*Node Output, nset=$(v)",
							"RF,"])
	end
	push!(stepString,"*End Step")
	return stepString
end

"""

	addStep!(abq::AbqModel, lc::LoadCase, out::Array{Output,1})

"""
function addStep!(abq::AbqModel,lc::LoadCase,out::Array{Output,1})
	push!(abq.steps,Step(lc.name,lc,out))
end

"""

	addStep!(abq::AbqModel, lc::LoadCase)

"""
function addStep!(abq::AbqModel,lc::LoadCase)
	push!(abq.steps,Step(lc.name,lc))
end
