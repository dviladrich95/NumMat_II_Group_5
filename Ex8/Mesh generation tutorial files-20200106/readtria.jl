"""
**Summary:**
```
xy,ele,vert,nvert,nele,nidv,ntop,nide=function readtriafile(fname::ASCIIString)
```

**Description:**

Reads a  *triangle* grid from the file `fname` and stores data in intermediate
structure consisting of

- `xyz::Array{Float64,2}` vertices positions
- `e2p::Array{Int64,2}`   element connectivity & IDs
- `idp::Array{Int64,2}`   vertex numbers & IDs
- `nvert::Int64`          number of vertices
- `nele::Int64`           number of elements
- `nidv::Int64`           number of vertex IDs
- `ntop::Int64`           number of entities in connectivity
- `nide::Int64`           number of element IDs

as read from *triangle* node and ele files.
"""

stringtoint(s)   = parse(Int64, s)
stringtofloat(s) = parse(Float64, s)

function readtria(fname::String)

  # read nodes
    fid = open(fname * ".node", "r")

    mesh_header = split(strip(readline(fid))) # read file header
    npoint      = stringtoint(mesh_header[1]) # number of points
    ndim        = stringtoint(mesh_header[2]) # dimension
    npattr      = stringtoint(mesh_header[3]) # number of point attributes
    nbound      = stringtoint(mesh_header[4]) # number of boundary markers

  # allocate
    x    = zeros(Float64, npoint, 1)
    y    = zeros(Float64, npoint, 1)
    idp = zeros(Int64, npoint, 1 + nbound)

  # set vertex positions and attributes
    for i = 1:npoint
        data = split(strip(readline(fid)))
        x[i] = stringtofloat(data[2])
        y[i] = stringtofloat(data[3])

        for j = 1:nbound
            idp[i,j] = stringtoint(data[3 + j])
        end
    end
    close(fid)

  # read elements
    fid = open(fname * ".ele", "r")
    mesh_header = split(strip(readline(fid))) # read header

    nele   = stringtoint(mesh_header[1]) # number of points
    nphi   = stringtoint(mesh_header[2]) # dimension
    neattr = stringtoint(mesh_header[3]) # number of element attributes

  # allocate
    e2p   = zeros(Int64, nele, nphi)
    ide   = zeros(Int64, nele, neattr)
  # set element connectivity and attributes
    for k = 1:nele
        data = split(strip(readline(fid)))
        for i = 1:nphi
            e2p[k,i] = stringtoint(data[i + 1])
        end
        for i = 1:neattr
            ide[k,i] = stringtoint(data[i + nphi + 1])
        end
    end
    close(fid)
    return x, y, npoint, nele, e2p, idp, ide
end