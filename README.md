# AbaqusUnitCell.jl
AbaqusUnitCell.jl provides functions for applying **periodic boundary conditions** on a **repeating unit cell** modelled in **Abaqus**.
The package contains functions for the following tasks:
- reading an Abaqus input-file,
- defining periodic boundary conditions,
- defining homogeneous displacement boundary conditions, and
- saving a new input-file with the defined boundary conditions.

**The following sections are currently in progress and curreently incomplete.**

## Model requirements
- The planes, which are defined by the periodic surfaces of the structure, should form a cuboid
- the surfaces have to be perpendicular or parallel to the coordinate axes
- Mesh has to be the same on opposing surfaces
- If the structure is cuboid, the vertices can be found automatically, otherwise, the vertices have to be modelled as Reference-Points in the Assembly
- When using tied surfaces or interactions:
	- if the coupled surfaces intersect with a periodic surface, the slave nodes have to be collected in a set called "Slaves"
	- These nodes will then be ignored in the equations to prevent overconstraints.

## Global and local coordinate system
For onedimensionally and twodimensionally periodic structures the applied boundary conditions depend on the directions of the periodicity.
While the structure is modelled in the **global** x-y-z-coordinate system, the periodic boundary conditions as well as the load boundary conditions are defined in a **local** I-II-III-coordinate system.
This distinction between global and local coordinates simplifies the definition of the PBCs.

Hereby, the local coordinate system is defined by choosing on axis of the global coordinate system as *reference axis*.
Direction I is then oriented in the same direction as the reference axis.
Similar to the global coordinates, the axis I, II, and III form an orthogonal right-hand-sided coordinate system.

| ref. axis | I   | II  | III |
| --------- | --- | --- | --- |
| x         | x   | y   | z   |
| y         | z   | x   | y   |
| z         | y   | z   | x   |

## Periodicity
In `AbaqusUnitCell.jl` onedimensional, twodimensional, or threedimensional periodicity can be defined for the 3D structure.
The global direction in which the structure is periodic is defined by the choice of the reference axis as follows:
- onedimensional periodicity: the structure is periodic in direction of the reference axis, which means, the coupled surfaces are perpendicular to the local direction I.
- twodimensional periodicity: the reference axis is perpendicular to the free surface. Therefore, the periodic surfaces are perpendicular to directions II and III.
- threedimensional periodicity: Each surface is couppled to the opposite one.

## Name convention for surfaces, edges and vertices
For the definition of the periodic boundary conditions the six surfaces of the repeating unit cell are named *North*, *South*, *East*, *West*, *Top*, and *Bottom*.
The assignment of the surface names depends on the definition of the local coordinate system.
*Top* and *Bottom* are always perpendicular to axis I, *East* and *West* are perpendicular to axis II, *North* and *South* are perpendicular to axis III (see figure ??).

Figure with name convention

The edges and vertices are named after the surfaces intersecting in this edge or vertex; for a shorter notation, only the first letters are taken into account (e.g. *SW* for the edge of *South* and *West* or **NWB** for the vertex of *North*, *West*, and *Bottom*).

## Periodic Boundary Conditions

## Load boundary conditions
The loads can be defined by the macroscopic strain tensor $\langle\varepsilon\rangle$.
This macroscopic strain tensor has to be formulated in the **local** coordinates!
For 1D periodicity:
$$\langle\varepsilon\rangle_{I\;I}$$
For 2D periodicity:
$$\langle\varepsilon\rangle = \begin{pmatrix}
\langle\varepsilon_{II\;II}\rangle & \langle\varepsilon_{II\;III}\rangle \\
\langle\varepsilon_{III\;II}\rangle & \langle\varepsilon_{III\;III}\rangle \\
\end{pmatrix}$$
For 3D periodicity:
$$\langle\varepsilon\rangle = \begin{pmatrix}
\langle\varepsilon_{I\;I}\rangle & \langle\varepsilon_{I\;II}\rangle & \langle\varepsilon_{I\;III}\rangle \\
\langle\varepsilon_{II\;I}\rangle & \langle\varepsilon_{II\;II}\rangle & \langle\varepsilon_{II\;III}\rangle \\
\langle\varepsilon_{III\;I}\rangle & \langle\varepsilon_{III\;II}\rangle & \langle\varepsilon_{III\;III}\rangle \\
\end{pmatrix}$$

It should be noted that only macroscopic strains in the direction of periodicity can be reasonably applied on the repeating unit cell. 

## Exported functions
### `AbqModel(path::String)`

### `setRefAxis!(abq::AbqModel, axis::String)`

### `setTolerance!(abq::AbqModel, tol::Number)`

### `setExceptions!(abq::AbqModel, exc::Array{String,1})`

### `setVertexFinder!(abq::AbqModel, val::Bool)`

### `setPBCdim!(abq::AbqModel, dim::Float64)`

### `nodeDesignation!(abq::AbqModel)`

### `pbc!(abq::AbqModel)`

### `LoadCase(strain::Array{<:Number,2}, abq::AbqModel)`
